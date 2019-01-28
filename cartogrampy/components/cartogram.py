import numpy as np
from numpy.linalg import norm
from geopandas import GeoDataFrame
import pandas as pd 
from shapely.geometry import asShape, Polygon
from pysal.weights.Contiguity import Rook, W
from geopandas import GeoDataFrame
from numpy import mean, sum, pi, power, maximum, minimum, array, concatenate, lexsort, any, argsort, searchsorted, split, take
from shapely.affinity import scale, translate

def borders_from_dataframe(df, idVariable=None,geom_field = 'geometry'):

    """Returns a PySAL weights object in which weights are lengths of shared border.

    Args:
        df (geopandas.geodataframe.GeoDataFrame): Input GeoDataFrame.
        idVariance (optional): Field name of values uniquely identifying rows in geodataframe, defaults to None.
        geom_field (str, optional): Field name of geometry column in input GeoDataFrame, defaults to 'geometry'.

    Returns:
        PySAL W weights object.
    """

    # Function to compute neighbors and length of shared borders
    # Based on: https://github.com/volaya/processing_pysal/blob/master/ext-libs/pysal/contrib/shared_perimeter_weights.py

    rook = Rook.from_dataframe(df,idVariable=idVariable)
    polygons = df.set_index(idVariable)[geom_field].apply(lambda x: asShape(x)).to_dict()
    new_weights = {}
    for i in rook.neighbors:
        a = polygons[i]
        new_weights[i] = [a.intersection(polygons[j]).length for j in rook.neighbors[i]]
    return W(rook.neighbors,new_weights)

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

        # make copy of the geodataframe containing only the relevant fields
        if self.id_field:
            gdf = self.gdf[[self.value_field, self.self.id_field, self.geom_field]].copy()
        else:
            gdf = self.gdf[[self.value_field, self.geom_field]].copy()

        # calculate geometry positions based on
        if position.lower() in ['centroid']:
            gdf['cent'] = gdf[self.geom_field].centroid
        elif position.lower() in ['center', 'centre']:
            gdf['cent'] = gdf[self.geom_field].envelope.centroid
        elif position.lower() in ['representative point', 'rep']:
            gdf['cent'] = gdf[self.geom_field].representative_point()
        else:
           # if position parameter not recognised, default the centroid
           print('position parameter invalid, using centroid.')
           gdf['cent'] = gdf[self.geom_field].centroid

        # work out the value densities and ranks
        gdf['density'] = gdf[self.value_field]/gdf.area
        gdf['rank'] = gdf['density'].rank(axis=0, method='first', ascending=False)

        # get appropriate anchor depending on whether anchor_id or anchor_rank  have been specified
        if anchor_id:
            if self.id_field:
                if anchord_id in gdf[self.id_field].values:
                    anchor = gdf[gdf[self.id_field] == anchor_id]['density'].values[0]
                else:
                    print("anchor_id not recognised in self.id_field, defaulting to anchor_rank = 1")
                    anchor = gdf[gdf['rank'] == 1]['density'].values[0]
            else:
                print("self.id_field not specified, for anchord_id, defaulting to anchor_rank = 1")
                anchor = gdf[gdf['rank'] == 1]['density'].values[0]
        else:
            anchor = gdf[gdf['rank'] == anchor_rank]['density'].values[0]

        # work out the scaling for each polygon
        gdf['scale'] = (1.0/power(anchor, 0.5)) * power(gdf[self.value_field]/gdf.area,0.5)

        # NB affine transformations are linear
        new_geom = [scale(g[1][self.geom_field], xfact=g[1]['scale'], yfact=g[1]['scale'],origin=g[1]['cent']) for g in gdf.iterrows()]

        # clean up
        del gdf['density'], gdf['rank'], gdf['cent']

        return GeoDataFrame(gdf, geometry=new_geom)

    def dorling(self,
              position='centroid',
              ratio=0.4,
              friction=0.25,
              iterations=99,
              verbose=True):

        if self.id_field:
            if position in self.gdf.columns.values:
                df = self.gdf[[self.id_field,
                               self.value_field,
                               self.geom_field,
                               position]].copy()
                id_field = self.id_field
            else:
                df = self.gdf[[self.id_field,
                               self.value_field,
                               self.geom_field,
                               position]].copy()
                id_field = self.id_field
        else:
            if position in self.gdf.columns.values:
                df = self.gdf[[self.value_field,
                               self.geom_field,
                               position]].copy()
                df['id_field'] = df.index
                id_field = 'id_field'
            else:
                df = self.gdf[[self.value_field,
                               self.geom_field]].copy()
                df['id_field'] = df.index
                id_field = 'id_field'

        # get lengths of shared borders
        wp = borders_from_dataframe(df, id_field, self.geom_field)

        # get dictionary of perimeters
        perimeter = df.set_index(id_field)[self.geom_field].length.to_dict()

        # now convert df to geodataframe of centers according to pos.
        if position in df.columns.values:
            try:
                df = GeoDataFrame(df, geometry=df[position])
            except:
                pass
        elif position.lower() == 'centroid':
            df = GeoDataFrame(df, geometry=df[self.geom_field].centroid)
        elif position.lower() in ['center', 'centre']:
            df = GeoDataFrame(df, geometry=df[self.geom_field].envelope.centroid)
        elif position.lower() in ['representative', 'rep']:
            df = GeoDataFrame(df, geometry=df[self.geom_field].representative_point())
        else:
            print("Did not recognize position argument, using centroid")
            df = GeoDataFrame(df, geometry=df[self.geom_field].centroid)


        # work out scale and radii - seems inefficient to deal with the geometries in this way.
        # no idea how to write this is a PEP8 compliant way.
        total_dist = sum([df.loc[df[id_field] == i,self.geom_field].centroid.tolist()[0].distance(df.loc[df[id_field] == j,self.geom_field].centroid.tolist()[0]) for i in wp.neighbors for j in wp[i]])
        total_radius = sum([power(df.loc[df[id_field] == i,self.value_field].values[0]/pi,0.5) + power(df.loc[df[id_field] == j,self.value_field].values[0]/pi,0.5) for i in wp.neighbors for j in wp[i]])

        scale = total_dist / total_radius

        radius = df.set_index(id_field)[self.value_field].apply(lambda x: scale * power(x/pi, 0.5)).to_dict()
        widest = radius[max(radius, key=lambda i: radius[i])]

        if verbose:
            print('Scaling by ', scale, ' widest is ', widest)

        # start the main loop
        for i in range(iterations):
            displacement = 0.0
            for idx, row in df.iterrows():
                # set up for values for each iteration
                xrepel = 0.0
                yrepel = 0.0
                xattract = 0.0
                yattract = 0.0
                closest = widest

                # get neighboring circles
                distband = row[self.geom_field].centroid.buffer(widest + radius[row[id_field]])
                neighbors = df[df[self.geom_field].centroid.within(distband)]
                # remove self from neighbors dataframe
                neighbors = neighbors[neighbors[id_field] != row[id_field]]

                # for neighbors calculate the attractive and repulsive forces acting
                if len(neighbors) > 0:

                    for nidx, nrow in neighbors.iterrows():
                        dist = row[self.geom_field].centroid.distance(nrow[self.geom_field].centroid)
                        if dist < closest:
                            closest = dist
                        overlap = radius[row[id_field]] + radius[nrow[id_field]] - dist
                        # calculate repelling forces for instersecting circles
                        if overlap > 0.0:
                            if dist > 1.0:
                                xrepel -= overlap * (nrow[self.geom_field].centroid.x - row[self.geom_field].centroid.x) / dist
                                yrepel -= overlap * (nrow[self.geom_field].centroid.y - row[self.geom_field].centroid.y) / dist

                        # calculate attractive forces for circles
                        if overlap < 0.0:
                            try:
                                overlap = abs(overlap) * wp[row[id_field]][nrow[id_field]]/perimeter[row[id_field]]
                            except KeyError:
                                overlap = 0.0
                            xattract = xattract + overlap * (nrow[self.geom_field].centroid.x - row[self.geom_field].centroid.x) / dist
                            yattract = yattract + overlap * (nrow[self.geom_field].centroid.y - row[self.geom_field].centroid.y) / dist

                # Calculate combined effect of attraction and repulsion
                atrdst = power(power(xattract, 2) + power(yattract, 2), 0.5)
                repdst = power(power(xrepel, 2) + power(yrepel, 2), 0.5)
                if repdst > closest:
                    xrepel = closest * xrepel / (repdst + 1.0)
                    yrepel = closest * yrepel / (repdst + 1.0)
                    repdst = closest
                if repdst > 0:
                    xtotal = (1.0 - ratio) * xrepel + ratio * (repdst * xattract / (atrdst + 1.0))
                    ytotal = (1.0 - ratio) * yrepel + ratio * (repdst * yattract / (atrdst + 1.0))
                else:
                    if atrdst > closest:
                        xattract = closest * xattract / (atrdst + 1.0)
                        yattract = closest * yattract / (atrdst + 1.0)
                    xtotal = xattract
                    ytotal = yattract
                displacement += power(power(xtotal, 2) + power(ytotal, 2), 0.5)

                # Record the vectors
                xvector = friction * xtotal
                yvector = friction * ytotal

                # update position
                df.loc[idx, self.geom_field] = translate(row[self.geom_field], xoff=xvector, yoff=yvector)

            displacement = displacement/len(df)

            if verbose:
                print("iter: ",i," displacement: ", displacement)

        df = df.merge(pd.DataFrame(radius.items(), columns=[id_field, 'Radius']),on = id_field)
        return GeoDataFrame(df,geometry=[df.loc[b,'geometry'].buffer(df.loc[b,'Radius']) for b in range(len(df))])

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
                print("iteration:", i+1, "Mean Error:", mean(maximum(gdf[self.geom_field].area,Desired)/minimum(gdf[self.geom_field].area,desired_area)),"Max Error:",max(maximum(gdf[self.geom_field].area,desired_area)/minimum(gdf[self.geom_field].area,desired_area)))

        return gdf
