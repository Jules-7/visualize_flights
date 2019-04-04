""" Flight visualization for inspection.

    Visualize a flown trajectory of a flight with its planned
    route with respect to MUAC AOR boundary.

    Table to query, as well as
    flight ifplid and date must be provided as an input.
"""

import pandas as pd
import matplotlib.pyplot as plt

from connect_to_db import db_connection
from set_map import basemap_setup


if __name__ == "__main__":

    m = basemap_setup()

    # AOR boundary
    aor = pd.read_csv('aor')
    x_muac, y_muac = m(aor['lon'].values, aor['lat'].values)
    plt.plot(x_muac, y_muac, color='royalblue', linewidth=2)

    conn = db_connection()
    cursor = conn.cursor()

    table_db = raw_input('Enter table to query:')
    ifpl_id = raw_input('Enter ifplid: ')
    target_date = raw_input('Enter date in format YYYY-MM-DD (for sbx2 only July 2018 is available): ')

    if 'atm' in table_db:
        flown_trajectory = pd.read_sql(""" select t.x lon, t.y lat, t.z alt, TO_DATE('19700101','yyyymmdd') + (t.w/24/60/60) t 
                            from (select * from (
                                            select * from %s where ifpl_id = '%s' and target_date = to_date('%s', 'YYYY-MM-DD') 
                                                    order by snapshot_time desc) 
                                            where rownum=1) f, table(sdo_util.getvertices(f.trajectory)) t""" % (table_db, ifpl_id, target_date), conn)
    else:
        flown_trajectory = pd.read_sql("""select t.x lon, t.y lat, t.z*100 alt, TO_DATE('19700101','yyyymmdd') + (t.w/1000/24/60/60) t
                                         from %s f, table(sdo_util.getvertices(f.trajectory)) t where ifplid = '%s'
                                         and target_date = to_date('%s', 'YYYY-MM-DD')""" % (table_db, ifpl_id, target_date), conn)

    if 'sbx2' in table_db:
        route = pd.read_sql("""select t.x lon, t.y lat from %s f, table(sdo_util.getvertices(f.route)) t
                                 where ifplid = '%s' and target_date = to_date('%s', 'YYYY-MM-DD')
                            """ % (table_db, ifpl_id, target_date), conn)

        # plot route
        route_latitudes = route['LAT'].values
        route_longitudes = route['LON'].values
        x, y = m(route_longitudes, route_latitudes)
        plt.plot(x, y, color='black', marker='o')

    conn.close()

    # plot flown trajectory
    flight_latitudes = flown_trajectory['LAT'].values
    flight_longitudes = flown_trajectory['LON'].values
    x, y = m(flight_longitudes, flight_latitudes)
    plt.plot(x, y, color='green')
    plt.plot(x[0], y[0], color='red', marker='o')

    plt.show()
