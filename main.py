import geopandas as gpd
import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from cartogrampy.components.cartogram import Cartogram

ldn_boro = gpd.read_file('data/LDN_Boro.geojson').to_crs(epsg=27700)
ldn_pop  = pd.read_excel('data/Pandas_Lon_Pop.xlsx')
data = ldn_boro.merge(ldn_pop, how='left', left_on='GSS_CODE', right_on='New Code')
c = Cartogram(data, value_field=2015 )
cc = c.dorling(iterations=200,verbose=True)
# Map the outcome
fig, ax = plt.subplots(1,2,figsize=(16,8))
ldn_boro.plot(column = 2015,cmap = 'PuRd', scheme = 'Quantiles', k=4,legend=True, ax=ax[0])
cc.plot(color='r',ax=ax[1])
ax[0].axis('equal')
ax[1].axis('equal')
