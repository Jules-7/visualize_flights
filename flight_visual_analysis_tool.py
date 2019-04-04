""" Flight Visual Analysis Tool

    This module provides with a GUI for flights visualisation.
    The GUI includes a set of filters for flights' selection and
    a map, where flights are displayed one-by-one. Additionally,
    the first aircraft heading and a list of bearings from the
    first position to each way point of the route are computed.

    Flights can be filtered by:
    - date
    - airdrome of destination (ADES)
    - airdrome of departure (ADEP)
    - upstream sector (NCP)
    - way point in the route (WPT)
    - entry type (lateral, vertical)
    - exit type (lateral, vertical)

    The map displays:
    - flown track
    - EFD route (route according to the flight plan with updates)
    - point of track correlation (i.e. the first point of the track)
    - point when flight is under control (by air traffic controller)

    The resultant plots can be saved as a figure and
    are meant for visual analysis only.
"""

import pandas as pd
import numpy as np
import Tkinter as tk
import matplotlib
import tkcalendar
import tkMessageBox
import ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from mpl_toolkits.basemap import Basemap
from connect_to_db import db_connection

matplotlib.use("TkAgg")


class DBConnection(object):
    """ Data Base connection """

    def __init__(self):
        self.db_connection_set()

    def db_connection_set(self):
        self.conn = db_connection()
        self.cursor = self.conn.cursor()

    def db_connection_close(self):
        self.conn.close()


class FlownTracks(object):
    """Query Data Base for flights
    satisfying user criteria"""

    def __init__(self, db=None, filters=None):
        self.db = db
        self.target_date = filters["target_date"]
        self.ncp = filters["ncp"]
        self.ades = filters["ades"]
        self.adep = filters["adep"]
        self.wpt = filters["wpt"]
        self.entry_type = filters["entry_type"]
        self.exit_type = filters["exit_type"]

    def flown_tracks_extract(self):
        """ Query all flights satisfying the conditions.
            Append conditions if provided."""

        flights_query = """select ifplid, target_date from flights 
                                        where target_date = to_date('{}', 'YYYY-MM-DD')
                                        and route_wpt is not null
                                        and ades is not null
                                        and adep is not null
                                        """.format(self.target_date)
        if self.ades:
            flights_query += """and ades = '{}'""".format(self.ades)

        if self.adep:
            flights_query += """and adep = '{}'""".format(self.adep)

        if self.ncp:
            flights_query += """and ncp = '{}'""".format(self.ncp)

        if self.wpt:
            flights_query += """and route_wpt like '%{}%'""".format(self.wpt)

        if self.entry_type:
            if self.entry_type == 'Lateral':
                # flights whose correlation point is >= FL245
                flights_query += """and corr_z >= 24500 * 0.3048"""
            elif self.entry_type == 'Vertical':
                # flights whose correlation point is < FL245
                flights_query += """and corr_z < 24500 * 0.3048"""

        if self.exit_type:
            if self.exit_type == 'Lateral':
                # flights whose correlation point is >= FL245
                flights_query += """and exit_z >= 24500 * 0.3048"""
            elif self.exit_type == 'Vertical':
                # flights whose correlation point is < FL245
                flights_query += """and exit_z < 24500 * 0.3048"""

        return self.db.cursor.execute(flights_query).fetchall()


class FlownTrack(object):
    """ """

    def __init__(self, db=None, target_date=None, ifplid=None, filters=None):
        self.db = db
        self.target_date = target_date
        self.ifplid = ifplid
        self.filters = filters

        self.EARTH_RADIUS = 6371e3  # [m], mean radius
        # ToDo discard route points that are further than X nm
        # from the track correlation point, (i.e. points in
        # Atlantic ocean, US, Eastern Europe, Asia, Africa)
        # self.MAX_DISTANCE_OF_INTEREST = None

    def trajectory_extract(self):
        flown_trajectory_query = """select t.x lon, t.y lat, t.z*100 alt, 
                                        TO_DATE('19700101','yyyymmdd') + (t.w/1000/24/60/60) t,
                                        f.under_control_x, f.under_control_y, ncp, ades, adep
                                        from flights f, table(sdo_util.getvertices(f.trajectory)) t
                                        where ifplid = '%s'
                                        and target_date = to_date('%s', 'YYYY-MM-DD')""" % (self.ifplid, self.target_date)
        self.trajectory = pd.read_sql(flown_trajectory_query, self.db.conn)

    def route_extract(self):
        """ For lateral entry and lateral exit use only the route way points, since only these
            are available for an upstream sector to give a direct-to.
            For vertical entries and exits and when conditions are not specified
            show also GEO points, which are provided in EFD messages by NM.
            This information might be useful for the climbing traffic. """
        if self.filters["entry_type"] == 'Lateral' and self.filters["exit_type"] == 'Lateral':
            # route without GEO points (e.g. only way points)
            route_query = """select t.x lon, t.y lat, f.route_wpt
                                     from flights f, table(sdo_util.getvertices(f.route)) t
                                     where ifplid = '%s'
                                     and target_date = to_date('%s', 'YYYY-MM-DD')""" % (self.ifplid, self.target_date)
        else:
            # route with GEO points (GEO points come in EFD messages from NM)
            route_query = """select t.x lon, t.y lat, f.route_wpt_with_geo route_wpt
                                     from flights f, table(sdo_util.getvertices(f.route_with_geo)) t
                                     where ifplid = '%s'
                                     and target_date = to_date('%s', 'YYYY-MM-DD')""" % (self.ifplid, self.target_date)
        route = pd.read_sql(route_query, self.db.conn)

        route_latitudes = route['LAT'].values
        route_longitudes = route['LON'].values
        route["WPT"] = route.loc[0, 'ROUTE_WPT'].split(',')

        flight_latitudes = self.trajectory['LAT'].values
        flight_longitudes = self.trajectory['LON'].values
        corr_lat = flight_latitudes[0]
        corr_lon = flight_longitudes[0]

        closest_segment_id = self.find_closest_segment(corr_lat, corr_lon, route_latitudes, route_longitudes)
        self.route = route.loc[closest_segment_id:].reset_index()

    def find_closest_segment(self, corr_lat, corr_lon, route_lats, route_longs):
        """ Finds the closest segment on the route to the correlation point

        :param corr_lat: latitude of the correlation point [DD]
        :param corr_lon: longitude of the correlation point [DD]
        :param route_lats: latitudes of the route points [DD]
        :param route_longs: longitudes of the route points [DD]
        :return: id of the route segment closest to the correlation point [-]
        """
        # ############## start test ###########################################
        # for test purposes compute intersection points and cross-distance for all points of the route
        # lat_1, lon_1 = route_lats[:-1], route_longs[:-1]
        # lat_2, lon_2 = route_lats[1:], route_longs[1:]

        # latitudes and longitudes of the closest point to p3 on each segment
        # closest_lats, closest_lons = self.nearest_point_on_segment(lat_1, lon_1, lat_2, lon_2, corr_lat, corr_lon)

        # compute distance to the closest point on each segment
        # dist_to_closest_points = self.distance(corr_lat, corr_lon, closest_lats, closest_lons)

        # returns same results as the statement above
        # compute cross-track distance from p3 to segment p1-p2
        # d_xt = self.distance_cross_segment(lat_1, lon_1, lat_2, lon_2, corr_lat, corr_lon)
        # print "cross dist to segment", d_xt.tolist()
        # ############## end test ###########################################

        # find the closest way point
        dist_to_wpts = self.distance(corr_lat, corr_lon, route_lats, route_longs)
        closest_wpt = np.argmin(dist_to_wpts)

        # select left and right segment
        # from the closest way point (i.e. the one before and the one after)
        wpts_ids = [closest_wpt - 1, closest_wpt, closest_wpt + 1]
        lat_1, lon_1 = route_lats[wpts_ids[:-1]], route_longs[wpts_ids[:-1]]
        lat_2, lon_2 = route_lats[wpts_ids[1:]], route_longs[wpts_ids[1:]]

        # latitudes and longitudes of the closest point to
        # the correlation point (p3) on each segment
        closest_lats, closest_lons = self.nearest_point_on_segment(lat_1, lon_1, lat_2, lon_2, corr_lat, corr_lon)

        # compute distance to the closest point on each segment
        dist_to_closest_points = self.distance(corr_lat, corr_lon, closest_lats, closest_lons)

        # d_xt results the same result as the statement above
        # compute cross-track distance from p3 to segment p1-p2
        # d_xt = self.distance_cross_segment(lat_1, lon_1, lat_2, lon_2, corr_lat, corr_lon)

        # check if p3 is on the segment
        on_segment = self.point_on_segment_check(lat_1, lon_1, lat_2, lon_2, corr_lat, corr_lon)

        if on_segment.size != 0:
            # return the first segment
            selected_segment_id = on_segment[0]

        else:
            # return the segment with the smallest distance
            selected_segment_id = np.argmin(dist_to_closest_points).item()

        closest_segment_id = wpts_ids[selected_segment_id]

        return closest_segment_id

    def point_on_segment_check(self, lat_1, lon_1, lat_2, lon_2, lat_3, lon_3):
        """ Check if a point lies on a segment

        :param lat_1: latitude of the segment start, p1, [DD]
        :param lon_1: longitude of the segment start, p1, [DD]
        :param lat_2: latitude of the segment end, p2, [DD]
        :param lon_2: longitude of the segment end, p2, [DD]
        :param lat_3: latitude of the point to be checked, p3, [DD]
        :param lon_3: longitude of the point to be checked, p3, [DD]
        :return: Boolean whether the point is on the segment or not
        """

        # compute distance from p1 to p2
        d12 = self.distance(lat_1, lon_1, lat_2, lon_2)

        # compute distance from p1 to p3
        d13 = self.distance(lat_1, lon_1, lat_3, lon_3)  # in meters

        # compute cross-track distance from p3 to segment p1-p2
        d_xt = self.distance_cross_segment(lat_1, lon_1, lat_2, lon_2, lat_3, lon_3)

        # d_xt < 1000 --> the cross-track distance is smaller than 1 kilometer
        # d13 < d12 --> along-track distance to the 'projected' point (with the cross-track distance)
        # is smaller then the length of the segment
        intersection = np.where((d_xt < 1000) & (d13 < d12))[0]

        return intersection

    def distance_cross_segment(self, lat_1, lon_1, lat_2, lon_2, lat_3, lon_3):
        """ Compute distance to a segment.
         Which is equivalent to the Cross-Track Distance (i.e. closest distance)
         Equations are taken from https://www.movable-type.co.uk/scripts/latlong.html

        :param lat_1: latitude of the segment start, p1, [DD]
        :param lon_1: longitude of the segment start, p1, [DD]
        :param lat_2: latitude of the segment end, p2, [DD]
        :param lon_2: longitude of the segment end, p2, [DD]
        :param lat_3: latitude of the point the distance is computed for, p3, [DD]
        :param lon_3: longitude of the point the distance is computed for, p3, [DD]
        :return: distance from p3 to segment p1p2 [m]
        """

        # compute distance from p1 to p3
        d13 = self.distance(lat_1, lon_1, lat_3, lon_3)  # in meters

        # compute bearing from p1 to p3
        b13, b13_deg = self.bearing(lat_1, lon_1, lat_3, lon_3)  # rad, deg

        # compute bearing from p1 to p2
        b12, b12_deg = self.bearing(lat_1, lon_1, lat_2, lon_2)

        # angular distance from p1 to p3
        d13_ang = d13 / self.EARTH_RADIUS

        # distance from p3 to segment p1-p2
        d_xt = np.arcsin(np.sin(d13_ang) * np.sin(b13 - b12)) * self.EARTH_RADIUS

        return np.abs(d_xt)

    def distance_along_segment(self, lat_1, lon_1, lat_2, lon_2, lat_3, lon_3):
        """ Distance along the segment from the start to the closest point
            on the segment to a third point
            (i.e. intersection with the cross-track distance)

        :param lat_1: latitude of the segment start, p1, [DD]
        :param lon_1: longitude of the segment start, p1, [DD]
        :param lat_2: latitude of the segment end, p2, [DD]
        :param lon_2: longitude of the segment end, p2, [DD]
        :param lat_3: latitude of the point the distance is computed for, p3, [DD]
        :param lon_3: longitude of the point the distance is computed for, p3, [DD]
        :return: distance along segment p1p2 to the intersection
                with the cross-track distance for point p3 [m]"""

        # compute distance from p1 to p3
        d13 = self.distance(lat_1, lon_1, lat_3, lon_3)  # in meters

        # angular distance from p1 to p3
        d13_ang = d13 / self.EARTH_RADIUS

        # compute cross-track distance from p3 to segment p1-p2
        d_xt = self.distance_cross_segment(lat_1, lon_1, lat_2, lon_2, lat_3, lon_3)

        # angular cross-track distance
        d_xt_ang = d_xt / self.EARTH_RADIUS

        d_at = np.arccos(np.cos(d13_ang) / np.cos(d_xt_ang)) * self.EARTH_RADIUS

        return d_at

    def nearest_point_on_segment(self, lat_1, lon_1, lat_2, lon_2, lat_3, lon_3):
        """Based on the cross-track distance and distance along
            the segment from the FIRST (i.e. START of the segment)
            point to the point of cross-track distance intersection,
            obtain the coordinates of the intersection
            (i.e. the nearest point on the segment)

        :param lat_1: latitude of the segment start, p1, [DD]
        :param lon_1: longitude of the segment start, p1, [DD]
        :param lat_2: latitude of the segment end, p2, [DD]
        :param lon_2: longitude of the segment end, p2, [DD]
        :param lat_3: latitude of the point the nearest point is computed for, p3, [DD]
        :param lon_3: longitude of the point the nearest point is computed for, p3, [DD]
        :return: latitude [DD] and longitude [DD]
                 of the nearest point to p3 on the segment p1p2 """

        # compute along-track distance to the closest point on the segment
        d_at = self.distance_along_segment(lat_1, lon_1, lat_2, lon_2, lat_3, lon_3)

        # compute segment bearing (from the first to the second point)
        b12, b12_deg = self.bearing(lat_1, lon_1, lat_2, lon_2)

        # estimate the intersection point from the distance and bearing starting at the first point
        lats, lons = self.point_from_distance_and_bearing(lat_1, lon_1, d_at, b12)

        return lats, lons

    def point_from_distance_and_bearing(self, lat_1, lon_1, dist, brng):
        """ Compute a point (coordinates) using a distance
            and a bearing to that point from the start point

        :param lat_1: latitude of the segment start, p1, [DD]
        :param lon_1: longitude of the segment start, p1, [DD]
        :param dist: distance to the point of interest from p1, [m]
        :param brng: bearing to the point of interest from p1, [rad]
        :return latitude [DD] and longitude [DD] of the point of interest
        """

        lat_1, lon_1 = np.radians(lat_1), np.radians(lon_1)

        # angular distance
        d_ang = dist / self.EARTH_RADIUS

        # compute destination point latitude
        lat_dest = np.arcsin(np.sin(lat_1) * np.cos(d_ang) +
                             np.cos(lat_1) * np.sin(d_ang) * np.cos(brng))

        # compute destination point longitude
        lon_dest = lon_1 + np.arctan2(np.sin(brng) * np.sin(d_ang) * np.cos(lat_1),
                                      np.cos(d_ang) - np.sin(lat_1) * np.sin(lat_dest))

        # convert from radians to degrees
        lat_dest = np.degrees(lat_dest)
        lon_dest = np.degrees(lon_dest)

        return lat_dest, lon_dest

    def distance(self, lat_1, lon_1, lat_2, lon_2):
        """ Great-circle distance computation using Haversine formula.
            Spherical Earth is assumed.

        :param lat_1: latitude of the first point, [DD]
        :param lon_1: longitude of the first point, [DD]
        :param lat_2: latitude of the second point, [DD]
        :param lon_2: longitude of the second point, [DD]
        :return distance between two points, [m]
        """
        lat_1, lon_1 = np.radians(lat_1), np.radians(lon_1)
        lat_2, lon_2 = np.radians(lat_2), np.radians(lon_2)
        delta_lat = lat_2 - lat_1
        delta_lon = lon_2 - lon_1

        a = np.sin(delta_lat / 2) * np.sin(delta_lat / 2) + \
            np.cos(lat_1) * np.cos(lat_2) * \
            np.sin(delta_lon / 2) * np.sin(delta_lon / 2)

        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

        d = self.EARTH_RADIUS * c

        return d

    def bearing(self, lat_1, lon_1, lat_2, lon_2):
        """ Compute bearing between two points.
            Bearing is computed in the range from -180 to +180.
            The conversion is performed to the range from 0 to 360
            (e.g. -90 degrees is converted into 270 degrees).

        :param lat_1: latitude of the first point, [DD]
        :param lon_1: longitude of the first point, [DD]
        :param lat_2: latitude of the second point, [DD]
        :param lon_2: longitude of the second point, [DD]
        :return bearing between two points, [rad], [deg]
        """

        lat_1, lon_1 = np.radians(lat_1), np.radians(lon_1)
        lat_2, lon_2 = np.radians(lat_2), np.radians(lon_2)
        delta_lon = lon_2 - lon_1

        y = np.sin(delta_lon) * np.cos(lat_2)
        x = np.cos(lat_1) * np.sin(lat_2) - np.sin(lat_1) * np.cos(lat_2) * np.cos(delta_lon)
        brng = np.arctan2(y, x)

        # if the angle in degrees is less than 0 (e.g. arctan2 returns values in -pi to +pi)
        # convert it to a positive angle
        brng_deg = np.degrees(brng)
        brng_deg[brng_deg < 0] = 360 + brng_deg[brng_deg < 0]
        brng = np.radians(brng_deg)

        return brng, brng_deg


class DisplayFlownTrack:
    """ GUI
    Main window includes widgets for flights' selection, map,
    flown trajectories and route bearing-related information.
    """
    def __init__(self, master):
        self.master = master

        self.flights_counter = 0  # used to keep track which flight to show next
        self.ERROR = False  # flag to raise if error occurs
        self.MIN_LON, self.MIN_LAT = -3, 46  # map lower left corner
        self.MAX_LON, self.MAX_LAT = 17, 57  # map upper right corner
        self.FEET_TO_M = 0.3048

        self.db = DBConnection()  # setup connection with the database
        self.load_available_ades()
        self.load_available_adep()
        self.load_available_ncp()
        self.load_available_wpt()

        self.create_widgets()   # layout and filling
        self.create_figure()
        self.events_handling()

    def events_handling(self):
        self.master.protocol("WM_DELETE_WINDOW", self.delete_window)  # capture delete window

    def delete_window(self):
        self.master.destroy()

    def load_available_ades(self):
        ades_query = """select distinct ades from flights where ades is not null"""
        self.ades_list = sorted([val[0] for val in self.db.cursor.execute(ades_query).fetchall()])

    def load_available_adep(self):
        adep_query = """select distinct adep from flights where adep is not null"""
        self.adep_list = sorted([val[0] for val in self.db.cursor.execute(adep_query).fetchall()])

    def load_available_ncp(self):
        ncp_query = """select distinct ncp from flights where ncp is not null"""
        self.ncp_list = sorted([val[0] for val in self.db.cursor.execute(ncp_query).fetchall()])

    def load_available_wpt(self):
        wpt_query = """select route_wpt from flights where route_wpt is not null"""
        wpt_result = [val[0] for val in self.db.cursor.execute(wpt_query).fetchall()]
        # remove first and last points (i.e. adep and ades)
        wpt_result_no_ad = [val.split(',')[1:-1] for val in wpt_result]
        self.wpt_list = sorted(list(set([item for sublist in wpt_result_no_ad for item in sublist])))

    def create_widgets(self):

        # ### Frames ##########################################################
        # the left one contains buttons, filters and displays table with numerical data
        # the right one contains map
        self.left_frame = tk.Frame(self.master, highlightbackground="red", highlightcolor="red", borderwidth=2)
        self.right_frame = tk.Frame(self.master)

        self.left_frame.pack(side='left', anchor='n')
        self.right_frame.pack(side='right')

        # ### Buttons #########################################################
        self.button_load = tk.Button(self.left_frame, text="Load flights", command=self.load_flights)
        self.button_next = tk.Button(self.left_frame, text="Next", command=self.show_next_flight)
        self.button_previous = tk.Button(self.left_frame, text="Previous", command=self.show_previous_flight)

        # ### Labels ##########################################################
        label_filters = tk.Label(self.left_frame, text="Filters", font=('Helvetica', 14, 'bold'))
        label_calendar = tk.Label(self.left_frame, text="Pick a date")
        label_calendar_note = tk.Label(self.left_frame, text="Data is available only for July 2018", fg='brown')
        label_ades = tk.Label(self.left_frame, text="ADES")
        label_adep = tk.Label(self.left_frame, text="ADEP")
        label_ncp = tk.Label(self.left_frame, text="NCP")
        label_wpt = tk.Label(self.left_frame, text="WPT")
        label_entry_type = tk.Label(self.left_frame, text="Entry type")
        label_exit_type = tk.Label(self.left_frame, text="Exit type")

        # ### Calendar ########################################################
        self.calendar = tkcalendar.DateEntry(self.left_frame, year=2018, month=7, day=1, width=8)

        # ### Comboboxes ######################################################
        ades_var = tk.StringVar(self.left_frame)
        adep_var = tk.StringVar(self.left_frame)
        ncp_var = tk.StringVar(self.left_frame)
        wpt_var = tk.StringVar(self.left_frame)
        ades_var.set(self.ades_list[0])
        ncp_var.set(self.ncp_list[0])
        wpt_var.set(self.wpt_list[0])

        self.combobox_ades = ttk.Combobox(self.left_frame, textvariable=ades_var, values=self.ades_list, width=8)
        self.combobox_adep = ttk.Combobox(self.left_frame, textvariable=adep_var, values=self.adep_list, width=8)
        self.combobox_ncp = ttk.Combobox(self.left_frame, textvariable=ncp_var, values=self.ncp_list, width=8)
        self.combobox_wpt = ttk.Combobox(self.left_frame, textvariable=wpt_var, values=self.wpt_list, width=8)

        # ### Radiobuttons ####################################################
        self.entry_var = tk.StringVar(self.left_frame)
        self.entry_var.set('Any')
        entry_options = ['Any', 'Lateral', 'Vertical']
        radiobutton_entry_a = tk.Radiobutton(self.left_frame, text=entry_options[0], variable=self.entry_var,
                                                value=entry_options[0])
        radiobutton_entry_l = tk.Radiobutton(self.left_frame, text=entry_options[1], variable=self.entry_var,
                                                value=entry_options[1])
        radiobutton_entry_v = tk.Radiobutton(self.left_frame, text=entry_options[2], variable=self.entry_var,
                                                value=entry_options[2])

        self.exit_var = tk.StringVar(self.left_frame)
        self.exit_var.set('Any')
        exit_options = ['Any', 'Lateral', 'Vertical']
        radiobutton_exit_a = tk.Radiobutton(self.left_frame, text=exit_options[0], variable=self.exit_var,
                                                value=exit_options[0])
        radiobutton_exit_l = tk.Radiobutton(self.left_frame, text=exit_options[1], variable=self.exit_var,
                                                value=exit_options[1])
        radiobutton_exit_v = tk.Radiobutton(self.left_frame, text=exit_options[2], variable=self.exit_var,
                                                value=exit_options[2])

        # ### Add widgets #####################################################
        label_filters.grid(row=1, column=1, sticky='w')

        label_calendar.grid(row=2, column=0)
        self.calendar.grid(row=2, column=1)
        label_calendar_note.grid(row=3, column=0, columnspan=2)

        label_ades.grid(row=4, column=0)
        label_adep.grid(row=5, column=0)
        label_ncp.grid(row=6, column=0)
        label_wpt.grid(row=7, column=0)

        self.button_load.grid(row=8, column=0, sticky='ew')
        self.button_next.grid(row=8, column=1, sticky='ew')
        self.button_previous.grid(row=8, column=2, sticky='ew')

        self.combobox_ades.grid(row=4, column=1)
        self.combobox_adep.grid(row=5, column=1)
        self.combobox_ncp.grid(row=6, column=1)
        self.combobox_wpt.grid(row=7, column=1)

        label_entry_type.grid(row=2, column=2)
        label_exit_type.grid(row=3, column=2)

        radiobutton_entry_a.grid(row=2, column=3)
        radiobutton_entry_l.grid(row=2, column=4)
        radiobutton_entry_v.grid(row=2, column=5)
        radiobutton_exit_a.grid(row=3, column=3)
        radiobutton_exit_l.grid(row=3, column=4)
        radiobutton_exit_v.grid(row=3, column=5)

    def create_figure(self):
        """Setup figure for basemap"""
        self.fig = Figure(figsize=(7, 7))
        self.ax = self.fig.add_subplot(111)

        # initialize the FigureCanvas to insert figure with basemap
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.right_frame)

        # toolbar is necessary for select, zoom, move and other buttons at the bottom
        toolbar = NavigationToolbar2Tk(self.canvas, self.right_frame)
        toolbar.update()

        self.canvas.get_tk_widget().pack(side=tk.BOTTOM)

        self.map_setup()

    def collect_filters(self):

        self.ERROR = False  # set error flag to false

        target_date = '20{2}-{0:0>2}-{1:0>2}'.format(*self.calendar.get().split('/'))

        # if target_date is outside the available range, raise a warning
        sel_year, sel_month = target_date.split('-')[:2]
        if int(sel_year) != 2018 or int(sel_month) != 7:
            tkMessageBox.showerror("Error", "Date is not available!\nChoose a date between 1 and 31 July 2018.")
            self.ERROR = True

        ades = self.combobox_ades.get().upper()
        adep = self.combobox_adep.get().upper()
        ncp = self.combobox_ncp.get().upper()
        wpt = self.combobox_wpt.get().upper()
        entry_type = self.entry_var.get()
        exit_type = self.exit_var.get()

        # for testing purposes
        # print 'filters'
        # print 'ades', ades
        # print 'adep', adep
        # print 'ncp', ncp
        # print 'wpt', wpt
        # print 'entry type', entry_type
        # print 'exit type', exit_type

        # collect all the filters
        self.filters = {'target_date': target_date,
                        'ades': ades,
                        'adep': adep,
                        'ncp': ncp,
                        'wpt': wpt,
                        'entry_type': entry_type,
                        'exit_type': exit_type}

    def load_flights(self):

        self.collect_filters()
        # reset the flights counter
        # (e.g. if new filters were applied and new flights are loaded)
        # set it to -1 so that with 'select_next_flight' method it becomes 0
        self.flights_counter = -1

        # if no errors occurred when collecting filters
        if not self.ERROR:

            flights = FlownTracks(self.db, self.filters)  # create flights object
            self.flights_ids = flights.flown_tracks_extract()  # get flights' ids

            # if no flights are found that satisfy the filtering criteria
            if not self.flights_ids:
                tkMessageBox.showerror("Error", "No flights satisfy selected criteria.\nSelect different filters.")
                return

            # display the first flight
            self.show_next_flight()

        return

    def map_setup(self):
        """Create Basemap and plot MUAC airspace boundaries"""
        llcrnrlon, llcrnrlat = self.MIN_LON, self.MIN_LAT  # lower left corner
        urcrnrlon, urcrnrlat = self.MAX_LON, self.MAX_LAT  # upper right corner
        lat_0 = (urcrnrlat - llcrnrlat) / 2.  # center latitude
        lon_0 = (urcrnrlon + abs(llcrnrlon)) / 2.  # center longitude

        # WGS-84: rsphere=(6378137.00, 6356752.3142)
        self.m = Basemap(projection='lcc', resolution='l',
                         llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                         rsphere=(6378137.00, 6356752.3142), lat_0=lat_0, lon_0=lon_0, ax=self.ax)
        self.m.fillcontinents(color='navajowhite', lake_color='lightsteelblue')
        self.m.drawcoastlines()
        self.m.drawcountries()
        self.m.drawmapboundary(fill_color='lightsteelblue')

        parallels = np.arange(llcrnrlat, urcrnrlat, 5.)
        # labels = [left,right,top,bottom]
        self.m.drawparallels(parallels, labels=[True, True, True, True])
        meridians = np.arange(llcrnrlon, urcrnrlon, 5.)
        self.m.drawmeridians(meridians, labels=[True, True, True, True])

        # Plot MUAC Airspace
        aor = pd.read_csv('aor')
        x_muac, y_muac = self.m(aor['lon'].values, aor['lat'].values)

        self.ax.plot(x_muac, y_muac, color='royalblue', linewidth=2)

    def show_next_flight(self):
        self.flights_counter += 1
        self.show_flight()

    def show_previous_flight(self):
        self.flights_counter -= 1
        self.show_flight()

    def show_flight(self):
        """Display a flight based on id"""
        # important to clear the previously displayed map
        self.ax.clear()

        # if available, extract id info on the next flight
        try:
            ifplid, target_date = self.flights_ids[self.flights_counter]

        except IndexError:  # if no more flights are available
            tkMessageBox.showerror("Warning", "No more flights satisfy selected criteria.\nSelect different filters.")
            return

        # convert datetime to a string for further queries
        target_date = target_date.strftime('%Y-%m-%d')

        # extract flight related information
        self.flight = FlownTrack(self.db, target_date, ifplid, self.filters)
        self.flight.trajectory_extract()
        self.flight.route_extract()

        # collect flight states based on altitude range (e.g. within or below MUAC)
        # take FL240 as the airspace bottom
        flight_within = self.flight.trajectory[self.flight.trajectory['ALT'] >= 24000 * self.FEET_TO_M]
        flight_below = self.flight.trajectory[self.flight.trajectory['ALT'] < 24000 * self.FEET_TO_M]

        # within
        flight_latitudes_within = flight_within['LAT'].values
        flight_longitudes_within = flight_within['LON'].values

        # below
        flight_latitudes_below = flight_below['LAT'].values
        flight_longitudes_below = flight_below['LON'].values

        # connect 'within' and 'below' arrays by adding a common point (if both exist)
        if flight_latitudes_below.size and flight_latitudes_within.size:

            below_first_index = self.flight.trajectory.index[self.flight.trajectory['LAT'] == flight_latitudes_below[0]]
            within_first_index = self.flight.trajectory.index[self.flight.trajectory['LAT'] == flight_latitudes_within[0]]

            below_last_index = self.flight.trajectory.index[self.flight.trajectory['LAT'] == flight_latitudes_below[-1]]
            within_last_index = self.flight.trajectory.index[self.flight.trajectory['LAT'] == flight_latitudes_within[-1]]

            # if 'below' part is first
            if below_first_index < within_first_index:
                # add the first point from the 'within' array
                flight_latitudes_below = np.append(flight_latitudes_below, [flight_latitudes_within[0]])
                flight_longitudes_below = np.append(flight_longitudes_below, [flight_longitudes_within[0]])
            else:
                # add the first point from the 'below' array
                flight_latitudes_within = np.append(flight_latitudes_within, [flight_latitudes_below[0]])
                flight_longitudes_within = np.append(flight_longitudes_within, [flight_longitudes_below[0]])

            # ToDo: flights with vertical entry and exit
            # it is possible, that aircraft climbed into MUAC, flew in MUAC for a while, and descended below MUAC
            # check this by looking if there are ids of 'below_index' that are larger then ids of 'within_index'
            # if that is the case, copy into a new array the ids that are larger than the last id in 'within'
            # and insert as the 0 index the last id in 'within'
            if below_last_index > within_last_index:
                pass

        # flight relevant information
        flight_times_ = self.flight.trajectory['T'].values
        correlation_time = pd.to_datetime(str(flight_times_[0])).strftime("%Y-%m-%d %H:%M:%S")
        correlation_lat = self.flight.trajectory['LAT'].values[0]
        correlation_lon = self.flight.trajectory['LON'].values[0]
        flight_ncp = self.flight.trajectory['NCP'].values[0]
        flight_ades = self.flight.trajectory['ADES'].values[0]
        flight_adep = self.flight.trajectory['ADEP'].values[0]

        # flight state at the moment of under control
        under_control_lat = self.flight.trajectory['UNDER_CONTROL_Y'].values
        under_control_lon = self.flight.trajectory['UNDER_CONTROL_X'].values

        # route coordinates
        route_latitudes = self.flight.route['LAT'].values
        route_longitudes = self.flight.route['LON'].values
        route_wpt_names = self.flight.route['WPT'].values

        # since self.ax has been cleaned, map needs to be re-created every time;
        # this results in a rather slow plots update
        self.map_setup()

        # convert latitudes and longitudes to map x,y values
        x_within, y_within = self.m(flight_longitudes_within, flight_latitudes_within)
        x_below, y_below = self.m(flight_longitudes_below, flight_latitudes_below)
        x_correlation, y_correlation = self.m(correlation_lon, correlation_lat)
        x_under_control, y_under_control = self.m(under_control_lon, under_control_lat)
        x_route, y_route = self.m(route_longitudes, route_latitudes)

        # plot route first so that it does not hide the other information
        self.ax.plot(x_route, y_route, color='black', marker='.',
                     markersize=12, label="route")  # route
        self.ax.plot(x_within, y_within, color='green', marker='.',
                     markersize=9, label="within MUAC")  # flown trajectory within MUAC
        self.ax.plot(x_below, y_below, color='orange', marker='.',
                     markersize=9, label="below MUAC")  # flown trajectory below MUAC
        self.ax.plot(x_correlation, y_correlation, color='red', marker='.',
                     markersize=10, label="correlation")  # track correlation point
        self.ax.plot(x_under_control[0], y_under_control[0], color='blue', marker='.',
                     markersize=10, label="assumed")  # under control point

        self.ax.legend(bbox_to_anchor=(1., 0.), ncol=3, loc='best',
                       markerscale=1, fontsize=9, borderaxespad=2)

        # plot way points' ids
        for i, wpt_name in enumerate(route_wpt_names):
            self.ax.annotate(i+1, xy=(x_route[i], y_route[i]))

        self.ax.set_title("IFPLID: %s. ADEP: %s. ADES: %s. \nNCP: %s. Correlation time %s\n" %
                          (ifplid, flight_adep, flight_ades, flight_ncp, correlation_time))

        # draw the map in the figure within the window frame
        self.canvas.draw()

        # compute the first bearing of the trajectory
        # and bearings from the correlation point to each way point of the route
        self.bearings_table()

    def bearings_table(self):
        """Plot a table containing the heading/bearing from
        the trajectory first point and bearings from
        the first point to all waypoints downstream."""
        bearing_start_row = 10  # starting from this row in the grid bearing info is displayed

        # clean previously displayed bearings
        for item in self.left_frame.grid_slaves():
            if int(item.grid_info()["row"]) >= bearing_start_row:
                item.grid_forget()

        # compute heading for each segment of the flown trajectory
        traj_points_count = self.flight.trajectory.index.size
        lat_1_t = self.flight.trajectory.loc[:traj_points_count-2, "LAT"].reset_index(drop=True)
        lon_1_t = self.flight.trajectory.loc[:traj_points_count-2, "LON"].reset_index(drop=True)
        lat_2_t = self.flight.trajectory.loc[1:, "LAT"].reset_index(drop=True)
        lon_2_t = self.flight.trajectory.loc[1:, "LON"].reset_index(drop=True)
        brng_true_t, brng_true_deg_t = self.flight.bearing(lat_1_t, lon_1_t, lat_2_t, lon_2_t)

        # compute heading from the correlation point to all way points (i.e. route points)
        route_points_count = self.flight.route.index.size
        lat_corr = self.flight.trajectory.loc[0, "LAT"]
        lon_corr = self.flight.trajectory.loc[0, "LON"]
        lat_1_r = pd.Series([lat_corr] * route_points_count)
        lon_1_r = pd.Series([lon_corr] * route_points_count)
        lat_2_r = self.flight.route["LAT"]
        lon_2_r = self.flight.route["LON"]
        brng_true_r, brng_true_deg_r = self.flight.bearing(lat_1_r, lon_1_r, lat_2_r, lon_2_r)

        route_wpt_names = self.flight.route['WPT'].values

        # ### Labels ################################################
        label_bearings_data = tk.Label(self.left_frame, text="Bearing Data", anchor='w', font=('Helvetica', 14, 'bold'))
        label_first_heading = tk.Label(self.left_frame, text='First heading', anchor='w', font=('Helvetica', 10, 'bold'))
        label_bearings = tk.Label(self.left_frame, text='Bearing to WPT', anchor='w', font=('Helvetica', 10, 'bold'))
        label_wpt_ids = tk.Label(self.left_frame, text='ID', anchor='w', font=('Helvetica', 10, 'bold'))
        label_wpt_names = tk.Label(self.left_frame, text='Name', anchor='w', width=7, font=('Helvetica', 10, 'bold'))

        # ### Entry #################################################
        entr_first_heading = tk.Entry(self.left_frame, width=10)
        entr_first_heading.insert(0, "%d" % brng_true_deg_t[0])

        # ### Add widgets ###########################################
        label_bearings_data.grid(row=bearing_start_row, column=1)
        label_first_heading.grid(row=bearing_start_row+1, column=0)
        label_bearings.grid(row=bearing_start_row+1, column=1)
        label_wpt_ids.grid(row=bearing_start_row+1, column=2)
        label_wpt_names.grid(row=bearing_start_row+1, column=3)

        entr_first_heading.grid(row=bearing_start_row+2, column=0)

        # traj_route_brng = list(zip(brng_true_deg_t[:brng_true_deg_r.size], brng_true_deg_r))

        # display the bearings
        for row, angle in enumerate(brng_true_deg_r):

            label_wpt_id = tk.Label(self.left_frame, text='%d'%(row+1), anchor='w', width=2)
            label_wpt_id.grid(row=row+bearing_start_row+2, column=2)

            label_wpt_name = tk.Label(self.left_frame, text='%s'%route_wpt_names[row], anchor='w', width=7)
            label_wpt_name.grid(row=row+bearing_start_row+2, column=3)

            entr_brng = tk.Entry(self.left_frame, width=10)
            entr_brng.insert(0, "%d" % angle)
            entr_brng.grid(row=row+bearing_start_row+2, column=1)


if __name__ == "__main__":
    root = tk.Tk()
    root.geometry("1250x650")
    root.title('Flights Visual Analysis Tool')
    app = DisplayFlownTrack(root)
    root.mainloop()
