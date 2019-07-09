import numpy as np
import pandas as pd
import geopandas as gpd

from scipy.spatial.distance import cdist
from pysal.lib.weights.contiguity import Rook, W
from shapely.geometry import Polygon
from shapely.affinity import scale, translate
from functools import partial

def borders_from_dataframe(gdf, idVariable = None, geom_field = 'geometry'):
    """Returns a pairwise pandas DataFrame of neighbouring zones with length of shared border as the weight.
    
    Args:
        gdf (gpd.GeoDataFrame): geoPandas geodataframe of zones.
        idVariable (optional): Column name of id variable, index used by default.
        geom_field (str, optional): Field name of geometry column in input gpd.GeoDataFrame, defaults to 'geometry'.
        
    Returns
        Pandas DataFrame
    """
    rook = Rook.from_dataframe(gdf, idVariable = idVariable)
    if idVariable:
        gdf = gdf.set_index(idVariable)
        weights = {idx: [gdf.loc[idx,geom_field].intersection(gdf.loc[nid,geom_field]).length 
                         for nid in neighbours] 
                   for idx, neighbours in rook.neighbors.items()}
    else:
        weights = {idx: [gdf.loc[idx,geom_field].intersection(gdf.loc[nid,geom_field]).length 
                         for nid in neighbours] 
                   for idx, neighbours in rook.neighbors.items()}
    return W(rook.neighbors,weights).to_adjlist()

def paired_distances(X, Y):
    """Pairwise distances for two arrays of coordinates
    
    Args:
        X: numpy array (n,2)
        Y: numpy array (n,2)
    
    Returns
        numpy array (n,1)
    """
    Z = X - Y
    norms = np.einsum('ij,ij->i', Z, Z)
    return np.sqrt(norms, norms)

def repel(x, row, xrepel, yrepel):
    if x['dist'] > 1.0:
        xrepel -= x['overlap'] * (x['geometry'].x - row['geometry'].x) / x['dist']
        yrepel -= x['overlap'] * (x['geometry'].y - row['geometry'].y) / x['dist']
    
    return (xrepel, yrepel)
        
def attract(x, idx, Wp, row, perimeter, xattract, yattract):
   
    if sum((Wp['focal'] == idx) & (Wp['neighbor'] == x.name)) == 1:
        x['overlap'] = abs(x['overlap']) * float(Wp[(Wp['focal'] == idx) & (Wp['neighbor'] == x.name)]['weight']) / perimeter[idx]
    
    xattract += x['overlap'] * (x['geometry'].x - row['geometry'].x) / x['dist']
    yattract += x['overlap'] * (x['geometry'].y - row['geometry'].y) / x['dist']
    
    return (xattract, yattract)

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

        # make copy of the geodataframe containing only the relevant fields

        geodf = self.gdf[[self.value_field, self.id_field, self.geom_field]].copy()

        # calculate geometry positions based on
        if position.lower() in ['centroid']:
            geodf['cent'] = geodf[self.geom_field].centroid
        elif position.lower() in ['center', 'centre']:
            geodf['cent'] = geodf[self.geom_field].envelope.centroid
        elif position.lower() in ['representative point', 'rep']:
            geodf['cent'] = geodf[self.geom_field].representative_point()
        else:
           # if position parameter not recognised, default the centroid
           print('position parameter invalid, using centroid.')
           geodf['cent'] = geodf[self.geom_field].centroid

        # work out the value densities and ranks
        geodf['density'] = geodf[self.value_field]/geodf.area
        geodf['rank'] = geodf['density'].rank(axis=0, method='first', ascending=False)

        # get appropriate anchor depending on whether anchor_id or anchor_rank  have been specified
        if anchor_id:
            if self.id_field:
                if anchor_id in geodf[self.id_field].values:
                    anchor = geodf[geodf[self.id_field] == anchor_id]['density'].values[0]
                else:
                    print("anchor_id not recognised in self.id_field, defaulting to anchor_rank = 1")
                    anchor = geodf[geodf['rank'] == 1]['density'].values[0]
            else:
                print("self.id_field not specified, for anchord_id, defaulting to anchor_rank = 1")
                anchor = geodf[geodf['rank'] == 1]['density'].values[0]
        else:
            anchor = geodf[geodf['rank'] == anchor_rank]['density'].values[0]

        # work out the scaling for each polygon
        geodf['scale'] = (1.0/np.power(anchor, 0.5)) * np.power(geodf[self.value_field]/geodf.area,0.5)

        # NB affine transformations are linear
        new_geom = [scale(g[1][self.geom_field], xfact=g[1]['scale'], yfact=g[1]['scale'],origin=g[1]['cent']) for g in geodf.iterrows()]

        # clean up
        del geodf['density'], geodf['rank'], geodf['cent']

        return gpd.GeoDataFrame(geodf, geometry=new_geom)

    def dorling(self,
              position='centroid',
              ratio=0.4,
              friction=0.25,
              iterations=99,
              verbose=True):

        # Get lengths of shared borders.
        Wp = borders_from_dataframe(self.gdf)
        perimeter = self.gdf.length
        
        # Get polygon centroids
        df = gpd.GeoDataFrame(self.gdf.drop(columns = self.geom_field), geometry = self.gdf.centroid)
        
        # Get total distance
        focal = (np.stack(Wp.merge(df[self.geom_field].map(np.array).to_frame(),
                                   left_on='focal',
                                   right_index=True).sort_index()[self.geom_field]))
        
        neighbour = (np.stack(Wp.merge(df[self.geom_field].map(np.array).to_frame(),
                                       left_on='neighbor',
                                       right_index=True).sort_index()[self.geom_field]))
        total_dist = np.sum(paired_distances(focal, neighbour))
    
        # Get total radius
        focal_rad = (Wp.merge(df[[self.value_field]], 
                              left_on = 'focal', 
                              right_index=True).sort_index()[self.value_field])
        neighbour_rad = (Wp.merge(df[[self.value_field]], 
                                  left_on = 'neighbor', 
                                  right_index=True).sort_index()[self.value_field])
        total_radius = np.sum((focal_rad / np.pi)**0.5 + (neighbour_rad / np.pi)**0.5)
        
        # Calculate scale
        scale = total_dist/total_radius
        
        # add radii
        df['radius'] = np.power(df[self.value_field]/np.pi,0.5) * scale
        WIDEST = df['radius'].max()
        
        # algorithm
        for i in range(iterations):
            displacement = 0.0
            # for each geometry in turn
            for idx, row in df.iterrows():
                xrepel = 0.0
                yrepel = 0.0
                xattract = 0.0
                yattract = 0.0
                closest = WIDEST
                
                # Get neighbours and calculate forces
                nb = df[df.distance(row[self.geom_field]).between(0, WIDEST + row['radius'], inclusive=False)].copy()
                if len(nb) > 0:
                    nb['dist'] = nb[self.geom_field].distance(row[self.geom_field])
                    closest = WIDEST if nb['dist'].min() > WIDEST else nb['dist'].min()
                    nb['overlap'] = (nb['radius'] + row['radius']) - nb['dist']
                    
                    for idy, rowy in nb.iterrows():
                        if rowy['overlap'] > 0.0:
                            xrepel, yrepel = repel(rowy, row, xrepel, yrepel)
                        else:
                            xattract, yattract = attract(rowy, idx, Wp, row, perimeter, xattract, yattract)
                # Calculate combined effect of attraction and repulsion
                atrdst = (xattract ** 2 + yattract ** 2) **0.5
                repdst = (xrepel ** 2 + yrepel ** 2) ** 0.5
                
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
                        
                displacement += (xtotal ** 2 + ytotal **2) ** 0.5
                
                # Record the vectors
                xvector = friction * xtotal
                yvector = friction * ytotal
                
                # update position
                df.loc[idx, self.geom_field] = translate(row[self.geom_field], xoff = xvector, yoff = yvector)
            
            displacement = displacement/len(df)
            if verbose:
                print(f"iter: {i}; displacement: {displacement:.2f}")
        
        return gpd.GeoDataFrame(data = df.drop(columns = ['geometry','radius']), geometry = df.apply(lambda x: x['geometry'].buffer(x['radius']), axis=1))
        
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