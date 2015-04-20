
"""

@author: af26

24 February 2015 --- Grid, Mpa and Larva classes separated out into 
separate files.

Uses velocity fields from the POLCOMS 1/6 x 1/9 degree model
to advect and diffuse larvae.

Written in Python 2.7 (though I don't think there is anything here which isn't
Python 3) because I couldn't find the basemap package for Python 3 on Windows.

Only tested on Windows using Anaconda python installation.

The data file path will need to be modified. Data file netcdf4.

Includes larval behaviour. Starting at bed, rising through water 
column, then sinking. 

POLCOMS netcdf arrays are indexed [k,j,i] ie [depth, lat, lon]


"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.path as mplPath
#from scipy.interpolate import interp1d
import shapefile
#from bngtolatlon import OSGB36toWGS84
import networkx as NX
from grid_class import Grid
from mpa_class import Mpa
from larva_class import Larva

run_dir = 'polcoms1991/Run_20150402/'

graph_output_dir = run_dir + 'Networkdata/'
track_output_dir = run_dir + 'Trackdata/'
mpa_name_file = open(run_dir + 'MPA_names.txt', 'r') 


nc_fileu = ('C:/Users/af26/PolcommModelData/' + 
            '1991/NOCL_S12run420_1991_UV.nc')
nc_filet = ('C:/Users/af26/PolcommModelData/' + 
            '1991/NOCL_S12run420_1991_TSz.nc')
nc_fidu = Dataset(nc_fileu, 'r')
nc_fidt = Dataset(nc_filet, 'r')

# larvae are released from MPA_SOURCE

# NUM_LARVAE are released at the start of each day for RELEASE_WINDOW days
# for a single release at time zero set RELEASE_WINDOW negative

NUM_LARVAE = 30
RELEASE_WINDOW = 30

STARTDAY = 32

SECONDS_IN_DAY = 60.0 * 60.0 * 24.0
RADIUS_OF_EARTH = 6378160.0
M_TO_DEGREE = 360.0 / (2.0 * np.pi * RADIUS_OF_EARTH)

DT = 3600.0
DTDAYS = DT / SECONDS_IN_DAY

KM = np.array([1.0, 1.0, 0.0002]) #constant diffusion coefficient m2/s

VERTICAL_INTERP = False
ANIMATE = False
SETTLING = True
DEATH = True

# constants for larval behaviour
# larvae released at bed, head upwards with increasing swim speeds up to 
# average age SWIMMAX. Then swimming gradually directed more downwards from 
# age DESCENDAGE up to DEADAGE.

SWIMSLOW = 0.0          #initial swimming speed
SWIMFAST = 0.003      # max swimming speed
SWIMSTART = 0.0         #age in days at which start swimming
SWIMMAX = 14.0          #average age in days at which max swimming speed is reached
DESCENDAGE = 21.0       # average age at which probability of heading down starts
                        # to increase
FULLDESCENDAGE = 42.0    # now fully heading down
MINSETTLEAGE = 30.0     # minimum age at which can settle given suitable 
                        # habitat
DEADAGE = 63.0          # Average age at which dead
# just set DEADAGE to a large value if larvae are not dying
if not DEATH:
    DEADAGE = 1000.0

if DEATH:
    if RELEASE_WINDOW > 0:
        NRUNDAYS = DEADAGE + RELEASE_WINDOW + 1
    else:
        NRUNDAYS = DEADAGE + 1
else:
    NRUNDAYS = 65
    
# bring constants together for passing to larva class

RUN_CONST = [SECONDS_IN_DAY, M_TO_DEGREE, DT, KM, VERTICAL_INTERP]
SWIM_CONST = [SWIMSLOW,SWIMFAST,SWIMSTART,SWIMMAX,DESCENDAGE,FULLDESCENDAGE,
              MINSETTLEAGE,DEADAGE]

def readVelocityData(nc_fid, n):
    '''
    netcdf Dataset -> 3 * masked array[i,j,k]
    Returns the 3-d  u and v velocity fields read from netcdf dataset 'nc_fid'.
    For the 'n'th day of the year. (Jan 1 = 1) 
    Vertical velocity w is set to zero - needs calculating from continuity.
    Velocities are stored as barotropic and baroclinic components. 
    Add for total (check this).
    
    '''

    u = nc_fidu.variables['U'][n-1,:,:,:]
    v = nc_fidu.variables['V'][n-1,:,:,:]
#    ub = nc_fidu.variables['UB'][n-1,:,:]
#    vb = nc_fidu.variables['VB'][n-1,:,:]
                    
    w = u * 0.0
    
    return u, v, w
    
def read_shapefile(filename):
    sf = shapefile.Reader(filename)
    shapes = sf.shapes()
    records = sf.records()
    return shapes, records

# helper functions
        
def group_settle(mpa_sprite_group, larva_object):
    settled = False
    for mpa in set(mpa_sprite_group):
        if mpa.settles(larva_object):
            settled = True
    return settled

def group_group_settle(larval_sprite_group, mpa_sprite_group):
    for larva in set(larval_sprite_group):
        if group_settle(mpa_sprite_group, larva):
            larval_sprite_group.remove(larva)
            settled_group.add(larva)  
            
def release_larvae(source, num, release_day):
    
    for mpa in mpa_group:
        
        if mpa.get_sitename() == source:
            nlarvae = 0
            bbox = mpa.get_bbox()
            points = mpa.get_points()
            path = mplPath.Path(points)
            while nlarvae < num:
                x = -200.0
                y = -200.0
                while not path.contains_point((x,y)):
                    x = np.random.uniform(bbox[0],bbox[2])
                    y = np.random.uniform(bbox[1],bbox[3])
            # check not on land
                i,j = gridt.get_index_ne(x,y)
                if not gridu.is_on_land(i,j):
                    larvae_group.add(Larva([x, y, -1.0], [0.0,0.0,0.0],
                                           source, release_day, gridt,
                                           RUN_CONST, SWIM_CONST))
                    nlarvae = nlarvae + 1    

def save_tracks_to_file(nc_outfile):

    nc_ofid = Dataset(nc_outfile, 'w')
    
    nl = len(larvae_dead) + len(settled_group) + len(larvae_outofarea)
    
    if DEATH:
        maxt = int((DEADAGE + 1) * SECONDS_IN_DAY / DT)
    else:
        maxt = int((NRUNDAYS) * SECONDS_IN_DAY / DT)
        
    time = nc_ofid.createDimension('time', maxt)
    nlarvae = nc_ofid.createDimension('nlarvae', nl)
    lon = nc_ofid.createVariable('longitude','f8',('nlarvae','time',))
    lat = nc_ofid.createVariable('latitude','f8',('nlarvae','time',))
    dep = nc_ofid.createVariable('depth','f8',('nlarvae','time',))
    rt = nc_ofid.createVariable('release day','i',('nlarvae',))
    fate = nc_ofid.createVariable('fate','S1',('nlarvae',))
    
    i = 0
    
    for larva in settled_group:
        x, y = larva.get_track()
        z = larva.get_depth_history()
        release_time = larva.get_release_day()
        state = 'S'
           
        lon[i,0:len(x)] = x[0:len(x)]
        lat[i,0:len(y)] = y[0:len(y)]
        dep[i,0:len(z)] = z[0:len(z)]
        rt[i] = release_time
        fate[i] = state
        
        i = i + 1
        
    for larva in larvae_dead:
        x, y = larva.get_track()
        z = larva.get_depth_history()
        release_time = larva.get_release_day()
        state = 'D'
           
        lon[i,0:len(x)] = x[0:len(x)]
        lat[i,0:len(y)] = y[0:len(y)]
        dep[i,0:len(z)] = z[0:len(z)]
        rt[i] = release_time
        fate[i] = state
        
        i = i + 1
    
    for larva in larvae_outofarea:
        x, y = larva.get_track()
        z = larva.get_depth_history()
        release_time = larva.get_release_day()
        state = 'L'
           
        lon[i,0:len(x)] = x[0:len(x)]
        lat[i,0:len(y)] = y[0:len(y)]
        dep[i,0:len(z)] = z[0:len(z)]
        rt[i] = release_time
        fate[i] = state
        
        i = i + 1
        
    nc_ofid.close()
    
    
G = NX.DiGraph()

    # read in and calculate the model grid variables
    
gridt = Grid(nc_fidt)
gridu = Grid(nc_fidu)
    
# loop over protected areas, releasing larvae from each
        
for line in mpa_name_file:
    MPA_SOURCE = line[0:-1]
    
    print MPA_SOURCE
    
    np.random.seed(1)
    
    # set up group of mpas
    
    mpa_group = set([])
    
    # offshore SAC
    shapes, records = read_shapefile('C:/Users/af26/Shapefiles/UK_SAC_MAR_GIS_20130821b/UK_SAC_MAR_GIS_20130821b/SCOTLAND_SAC_OFFSHORE_20121029_SIMPLE3')
    for i in range(len(shapes)):
        mpa_group.add(Mpa(shapes[i], records[i],'OFF_SAC'))
        
    # SAC with marine components
    shapes, records = read_shapefile('C:/Users/af26/Shapefiles/UK_SAC_MAR_GIS_20130821b/UK_SAC_MAR_GIS_20130821b/SCOTLAND_SACs_withMarineComponents_20130821_SIMPLE3')
    for i in range(len(shapes)):
        mpa_group.add(Mpa(shapes[i], records[i],'MAR_SAC'))
        
    # Nature conservation MPA
    shapes, records = read_shapefile('C:/Users/af26/Shapefiles/MPA_SCOTLAND_ESRI/MPA_SCOTLAND_SIMPLE3')
    for i in range(len(shapes)):
        mpa_group.add(Mpa(shapes[i], records[i],'MPA'))
        
    # Irish SACs
    shapes, records = read_shapefile('C:/Users/af26/Shapefiles/SAC_ITM_WGS84_2015_01/SAC_Offshore_WGS84_2015_01')
    for i in range(len(shapes)):
        mpa_group.add(Mpa(shapes[i], records[i],'IRISH'))
        
    # read in the opening day's data
    
    u, v, w = readVelocityData(nc_fidu,STARTDAY)
    
    # initialise larvae. 
    # Using grids of larvae at the same depth around a central point.
    
    larvae_group = set([])
    larvae_outofarea = set([])
    settled_group = set([])
    larvae_dead = set([])

    # seed larvae randomly in a particular mpa
    
    release_larvae(MPA_SOURCE,NUM_LARVAE, 0)

# the main program loop            
            
    day = STARTDAY
    nsteps = int(SECONDS_IN_DAY * NRUNDAYS / DT)
    endday = False
    runtime = 0.0
    for t in range(nsteps):
        
        rundays = runtime / SECONDS_IN_DAY
        
    #    print runtime, rundays
        
        for larva in set(larvae_group):
            left = larva.update(DT, rundays, gridu, gridt, u, v, w)
            if left:
                larvae_outofarea.add(larva)
                larvae_group.remove(larva)
            else:
                if larva.dead():
                    larvae_dead.add(larva)
                    larvae_group.remove(larva)
        runtime = (t + 1) * DT
        # read in a new current data field at the start of each day (daily mean fields)
        
        if ((runtime % SECONDS_IN_DAY) < DT/2.0):
            day = day + 1
#            print day
            u, v, w = readVelocityData(nc_fidu,day)
            
        # release a new batch of larvae if still in RELEASE_WINDOW
            if int(round(rundays)) < RELEASE_WINDOW:
                release_larvae(MPA_SOURCE, NUM_LARVAE, rundays)
    
#            if ANIMATE:
#                plot_animate()
                
        if (SETTLING and (rundays > MINSETTLEAGE)):
            group_group_settle(larvae_group, mpa_group)
            
    # output the connectivity graph
    
    
    # build graph
    G.clear()
    
    G.add_node(MPA_SOURCE)
    
    for mpa in mpa_group:
        nsettled = mpa.get_settled()
        if nsettled != 0:
            mpa_name = mpa.get_sitename()
            weight = float(nsettled) / float(NUM_LARVAE)
            G.add_weighted_edges_from([(MPA_SOURCE,mpa_name,weight)])
                       
    # output graph to file
#
    outfile = open(graph_output_dir + MPA_SOURCE + '.graphml', 'w')
    NX.write_graphml(G,outfile)
    outfile.close()
    
    # save tracks to file
    
    track_outfile = track_output_dir + MPA_SOURCE + '.nc'

    save_tracks_to_file(track_outfile)

    

