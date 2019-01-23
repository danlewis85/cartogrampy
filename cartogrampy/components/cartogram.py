from numpy import mean, pi, power, maximum, minimum, array, concatenate, lexsort, any, str, argsort, searchsorted, split, take
from numpy.linalg import norm
from shapely.geometry import Polygon
from geopandas import GeoDataFrame

class Cartogram:
    def __init__(self,
               gdf,
               value_field,
               id_field   = None,
               geom_field = 'geometry'):
       
        # dataframe and column descriptors
        self.gdf         = gdf
        self.value_field = value_field
        self.id_field    = id_field
        self.geom_field  = geom_field
        
        
    def noncont(self, 
              position='centroid',
              anchor_rank=1,
              anchor_id=None):
        pass
        
    def dorling(self,
              position='centroid',
              ratio=0.4,
              friction=0.25,
              iterations=99,
              verbose=True):
        pass
    
    def dcn(self,
          iterations=99,
          verbose=True):
        
        #+TODO: see if this can be generalised in a method
        #+TODO: copy instead of view
        # set out the dataframe
        gdf = self.gdf[[self.geom_field,
                        self.value_field]].copy()
        #+TODO: Note that we might as well keep the index insted of creating a new id field
        if self.id_field: 
            id_field = self.id_field
            gdf[id_field] = self.gdf[id_field]
        else:
            id_field = 'id_field'
            gdf[id_field] = gdf.index
        
        # compute sum of value_field
        total_val = gdf[self.value_field].sum()
       
        # main loop
        for i in range(iterations):
            # calculate current stats for gdf
            # total area
            total_area = gdf[self.geom_field].area.sum() 
            # desired areas relative to total area of current iteration
            desired_area = total_area * (gdf[self.value_field]/total_val)
            # circular radius for each geometry as a function of area of current iteration
            radius = power(gdf[self.geom_field].area/pi, 0.5)
            # mass as difference in desired and actuall radii
            mass = power(desired_area/pi, 0.5) - radius
            # ratio of desired to current area as an error rate 
            size_error = maximum(gdf[self.geom_field].area,
                                  desired_area) / minimum(
                                      gdf[self.geom_field].area, desired_area)
            # work out force reduction factor based on mean error 
