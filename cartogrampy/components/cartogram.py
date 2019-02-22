import numpy as np
import pandas as pd
import geopandas as gpd

from scipy.spatial.distance import cdist
#from pysal.weights.Contiguity import Rook, W
from shapely.geometry import asShape, Polygon
from shapely.affinity import scale, translate
from functools import partial

def borders_from_dataframe(df, idVariable=None,geom_field = 'geometry'):

    """Returns a PySAL weights object in which weights are lengths of shared border.

    Args:
        df (geopandas.geodataframe.gpd.GeoDataFrame): Input GeoDataFrame.
        idVariance (optional): Field name of values uniquely identifying rows in geodataframe, defaults to None.
        geom_field (str, optional): Field name of geometry column in input gpd.GeoDataFrame, defaults to 'geometry'.

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
               id_field = None,
               geom_field = 'geometry'):
       
        # dataframe and column descriptors
        self.gdf         = gdf
        self.value_field = value_field
        self.id_field    = id_field
       
        # TODO: make index of frame id column if one exists
        # TODO: use the index for any id needs -> not filling in new cols 
        # create an id_field
        if not id_field:
            self.gdf["id_field"] = self.gdf.index
            self.id_field = "id_field"
            
        self.geom_field  = geom_field

        self.multi = any(gdf.geom_type == "MultiPolygon")

    @classmethod
    def multi2single(self, gdf, geom_field = 'geometry'):

        """Returns a GeoDataFrame in which multipart features are broken into singlepart features.

        Args:
            gdf (geopandas.geodataframe.gpd.GeoDataFrame): Input GeoDataFrame.
            geom_field (str, optional): Field name of geometry column in input gpd.GeoDataFrame, defaults to 'geometry'.

        Returns
            Geopandas GeoDataFrame

        """
        # split rows into single and multipart features
        gdf_single = gdf[gdf[geom_field].type == 'Polygon']
        gdf_multi = gdf[gdf[geom_field].type == 'MultiPolygon']
        # Apply function to separate multipolygons
        partial_sep = partial(self.__separate, geom_field = 'geometry')
        sep = gdf_multi.apply(partial_sep, axis=1).tolist()
        sep.append(gdf_single)
        # Join all singlepart features together
        out = pd.concat(sep).reset_index(drop = True)
        # assign crs
        out.crs = gdf.crs
        return out

    def noncont(self,
              position='centroid',
              anchor_rank=1,
              anchor_id=None):
        # NOTE: BROKEN
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
        gdf['scale'] = (1.0/np.power(anchor, 0.5)) * np.power(gdf[self.value_field]/gdf.area,0.5)

        # NB affine transformations are linear
        new_geom = [scale(g[1][self.geom_field], xfact=g[1]['scale'], yfact=g[1]['scale'],origin=g[1]['cent']) for g in gdf.iterrows()]

        # clean up
        del gdf['density'], gdf['rank'], gdf['cent']

        return gpd.GeoDataFrame(gdf, geometry=new_geom)

    def dorling(self,
              position='centroid',
              ratio=0.4,
              friction=0.25,
              iterations=99,
              verbose=True):
        # NOTE: BROKEN
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
                df = gpd.GeoDataFrame(df, geometry=df[position])
            except:
                pass
        elif position.lower() == 'centroid':
            df = gpd.GeoDataFrame(df, geometry=df[self.geom_field].centroid)
        elif position.lower() in ['center', 'centre']:
            df = gpd.GeoDataFrame(df, geometry=df[self.geom_field].envelope.centroid)
        elif position.lower() in ['representative', 'rep']:
            df = gpd.GeoDataFrame(df, geometry=df[self.geom_field].representative_point())
        else:
            print("Did not recognize position argument, using centroid")
            df = gpd.GeoDataFrame(df, geometry=df[self.geom_field].centroid)


        # work out scale and radii - seems inefficient to deal with the geometries in this way.
        # no idea how to write this is a PEP8 compliant way.
        total_dist = np.sum([df.loc[df[id_field] == i,self.geom_field].centroid.tolist()[0].distance(df.loc[df[id_field] == j,self.geom_field].centroid.tolist()[0]) for i in wp.neighbors for j in wp[i]])
        total_radius = np.sum([np.power(df.loc[df[id_field] == i,self.value_field].values[0]/np.pi,0.5) + np.power(df.loc[df[id_field] == j,self.value_field].values[0]/np.pi,0.5) for i in wp.neighbors for j in wp[i]])

        scale = total_dist / total_radius

        radius = df.set_index(id_field)[self.value_field].apply(lambda x: scale * np.power(x/np.pi, 0.5)).to_dict()
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
                atrdst = np.power(np.power(xattract, 2) + np.power(yattract, 2), 0.5)
                repdst = np.power(np.power(xrepel, 2) + np.power(yrepel, 2), 0.5)
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
                displacement += np.power(np.power(xtotal, 2) + np.power(ytotal, 2), 0.5)

                # Record the vectors
                xvector = friction * xtotal
                yvector = friction * ytotal

                # update position
                df.loc[idx, self.geom_field] = translate(row[self.geom_field], xoff=xvector, yoff=yvector)

            displacement = displacement/len(df)

            if verbose:
                print("iter: ",i," displacement: ", displacement)

        df = df.merge(pd.DataFrame(radius.items(), columns=[id_field, 'Radius']),on = id_field)
        return gpd.GeoDataFrame(df,geometry=[df.loc[b,'geometry'].buffer(df.loc[b,'Radius']) for b in range(len(df))])

    def dcn(self,
          iterations=99,
          verbose=True):
        
        # If id field is specified get a copy of the geodataframe with just the id,
        # geometry and value fields.
        geodf = self.gdf[[self.geom_field, self.value_field, self.id_field]].copy()

        # compute sum of value_field and store
        totalValue = geodf[self.value_field].sum()

        def _calc_factors(frame):

            if self.multi:
                # Dissolve singlepart back to multipart
                frame = frame.dissolve(by=self.id_field, as_index=False)

            # Desired areas relative to total area of current iteration.
            desired = totalArea * (frame[self.value_field]/totalValue)
            # Circular radius for each geometry as a function of area of current iteration
            radius = (frame.area/np.pi).pow(0.5)
            # Mass as difference in desired and actual radii
            mass = (desired/np.pi).pow(0.5) - radius
            # Ratio of desired to current area as an error rate.
            sizeError = np.maximum(frame.area,desired)/np.minimum(frame.area,desired)
            # Work out force reduction factor based on mean error.
            forceReductionFactor = 1.0 / (1.0 + sizeError.mean())
            # Create centroids for convenience. Need to be multi-polygon centroids
            # if multipolygon.
            centroids = frame.centroid

            return radius, desired, mass, sizeError, forceReductionFactor, centroids


        # Main loop for iterations.
        for i in range(iterations):
            
            # Total area of all geometries for current iteration
            totalArea = geodf.area.sum()

            # TODO: make sure this is inside the function below.
            if self.multi:
                # Now prepare the geodataframes - ensure singlepart polygons and
                # deduplicate geometry points.
                geodf = self.multi2single(geodf, self.geom_field)

            (radius,
             desired,
             mass,
             sizeError,
             forceReductionFactor,
             centroids) = _calc_factors(geodf)

            # TODO: Sort this out in terms of dataframes if possible.
            # Now that we have singlepart features, deduplicate geometry points.
            # Get polygon coords as Series of numpy arrays
            pnts = geodf.exterior.map(np.array)
            # create a lookup based on count of points per polygon
            cnts = pnts.map(len)
            # Concatenate all points into a single long list of points.
            pnts = np.concatenate(pnts.values)
            # Get indices of sorted points
            ind = np.lexsort(pnts.T)
            # Get unique points by comparing sorted array of points
            upnts = pnts[ind[np.concatenate(([True],np.any(pnts[ind[1:]]!=pnts[ind[:-1]],axis=1)))]]

            # Make the arrays into structured arrays for index lookup creation
            a = pnts.ravel().view([('x',np.float64),('y',np.float64)])
            b = upnts.ravel().view([('x',np.float64),('y',np.float64)])
            sorter = np.argsort(b)
            idv = sorter[np.searchsorted(b, a, sorter=sorter)]

            centxy = centroids.map(np.array)
            upnts0 = upnts.copy()

            def _process(changed,
                         original,
                         c1,c2,c3,
                         mask,
                         cxy,
                         dim=0):
                """
                abstraction to process the different dimension within the
                next loop
                """
                diff = original[:, dim] - cxy[dim]

                calc = lambda const: diff * (const * c1) + changed[:, dim]
                pos, neg  = calc(c2), calc(c3)

                changed[:,dim] = np.where(mask, pos, neg)

                return changed

            for idx, cxy in enumerate(centxy):
                # make distance vector
                dist = cdist(upnts,cxy.reshape((1,2)))[:,0]

                # create boolean filter
                mask = dist > radius[idx]

                # maths
                # constants that don't need to be done in _process
                # but depend on the loop
                # real meat and potatoes of the algorithm
                c1 = forceReductionFactor / dist
                c2 = mass[idx] * radius[idx] / dist
                c_ = dist/radius[idx]
                c3 = mass[idx] * c_**2 * (4 - (3 * c_))

                # prefil a function
                upnts = _process(upnts, upnts0, c1, c2, c3, mask, cxy, dim=0)
                upnts = _process(upnts, upnts0, c1, c2, c3, mask, cxy, dim=1)

            # Reconstruct the full points list from the unique points using take and the index list.
            # Remake the array - split by counts lookup.
            # NB without the -1 index you get a 0 length array at the end.
            repnts = np.split(np.take(upnts,idv,axis=0),cnts.cumsum())[:-1]
            new_geom = [Polygon(p) for p in repnts]
            geodf = gpd.GeoDataFrame(geodf[[self.value_field,self.id_field]],geometry=new_geom)

            if self.multi:
                # Dissolve back to original multipart polygon.
                geodf = geodf.dissolve(by = self.id_field,as_index=False)

            # TODO: sort out this mess. good start though.
            # TODO: use f strings
            if verbose:
                mean_error = np.mean(np.maximum(geodf[self.geom_field].area,desired)/np.minimum(geodf[self.geom_field].area,desired))
                max_error = max(np.maximum(geodf[self.geom_field].area,desired)/np.minimum(geodf[self.geom_field].area,desired))
                print("iteration: {}; Mean Error: {:.3f}; Max Error: {:.3f}".format(i+1, mean_error, max_error))

        return geodf
    @staticmethod
    def __separate(row, geom_field):
        """Helper function for _multi2single.

        Args:
            row (pd.Series): row objects from dataframe
            geom_field (str, optional): Field name of geometry column in input gpd.GeoDataFrame, defaults to 'geometry'.

        Returns
            Pandas DataFrame
        """
        df = pd.concat([gpd.GeoDataFrame(row).T] * len(row[geom_field]),ignore_index = True)
        df[geom_field] = row[geom_field]
        return df
