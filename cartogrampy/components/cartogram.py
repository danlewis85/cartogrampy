from numpy import mean, pi, power, maximum, minimum, array, concatenate, lexsort, any, argsort, searchsorted, split, take
import numpy as np
from numpy.linalg import norm
from shapely.geometry import Polygon
from geopandas import GeoDataFrame
import pandas as pd 
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
            force_reduction = 1.0 / (1.0 + mean(size_error))
            # create centroids for convenience. need to be multi-polygon centroids
            centroids = gdf[self.geom_field].centroid
            
            # if multipolygons exist convert to polygons
            multi = False
            #+TODO: Might be worth just seperating out the fields. Keeping them in a df
            # rather than just column variables just adds text
            if gdf[gdf[self.geom_field].type == 'MultiPolygon'][self.geom_field].any():
                data = {col:[] for col in gdf.columns if col != self.geom_field}
               
                
                # now iterate over gdf to get data and geometry rows
                geom = []
                for idx,row in gdf.iterrows():
                    if row[self.geom_field].type == 'Polygon':
                        for col in data.keys():
                            data[col].append(row[col])
                        geom.append(row[self.geom_field])
                    elif row[self.geom_field].type == 'MultiPolygon':
                        for g in row[self.geom_field]:
                            for col in data.keys():
                                data[col].append(row[col])
                            geom.append(g)
                    else:
                        print("Geometry Error")
                    
                # reset gdf to the singlepart feature geodataframe
                gdf = GeoDataFrame(data,geometry = geom)        
                # flag that the geometry in question is originally multipart.
                multi = True
            
            # deduplicate geometry points
            # create a series of numpy arrays of boundary points
            pnts = gdf[self.geom_field].apply(lambda x: array(x.exterior))
          
            
            # create a lookup based on count of points per polygon
            cnts = pnts.apply(lambda x: len(x))

            # Concatenate all points into a single long list of points.
            pnts = concatenate(array(pnts))

            # Calculate unique points in full dataset.
            # From: http://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
            ind = lexsort(pnts.T)
            upnts = pnts[ind[concatenate(([True],any(pnts[ind[1:]]!=pnts[ind[:-1]],axis=1)))]]

            # This creates an id lookup for pnts to upnts.
            # It creates string views of the data, which is WAY quicker to deal with that 2d arrays.
            # Using http://stackoverflow.com/questions/32191029/getting-the-indices-of-several-elements-in-a-numpy-array-at-once
            a = pd.Series(pnts.ravel()).astype(str)
            b = pd.Series(upnts.ravel()).astype(str)
            sorter = argsort(b)
            idv = sorter[searchsorted(b, a, sorter=sorter)]

            # This is the main algorithm, and the bit that really takes the time.
            points = []
            for xy in upnts:
                x = x0 = xy[0]
                y = y0 = xy[1]
                idx = 0
                for cxy in centroids.apply(lambda x: [x.x,x.y]):
                    dist = norm(array(cxy)-array([x,y]))
                    if dist > radius[idx]:
                        fij = mass[idx] * radius[idx] / dist
                    else:
                        xf = dist / radius[idx]
                        fij = mass[idx] * (power(xf,2)) * (4 - (3 * xf))
                    fij = fij * force_reduction / dist

                    x = (x0 - cxy[0]) * fij + x
                    y = (y0 - cxy[1]) * fij + y
                    idx +=1
                points.append([x,y])
            points = array(points)
            # Reconstruct the full points list from the unique points using take and the index list.
            # Remake the array - split by counts lookup.
            # NB without the -1 index you get a 0 length array at the end.
            repnts = split(take(points,idv,axis=0),cnts.cumsum())[:-1]
            new_geom = [Polygon(p) for p in repnts]
            gdf = GeoDataFrame(gdf[[self.value_field,id_field]],geometry=new_geom)

            if multi:
                # Dissolve back to original multipart polygon.
                gdf = gdf.dissolve(by = id_field,as_index=False)
            if verbose:
                print("iteration:", i+1, "Mean Error:", mean(maximum(gdf[self.geom_field].area,Desired)/minimum(gdf[self.geom_field].area,desired_area)),"Max Error:",max(maximum(gdf[geom_field].area,desired_area)/minimum(gdf[self.geom_field].area,desired_area)))

        return gdf
