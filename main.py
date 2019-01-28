import geopandas as gpd
import pandas as pd
import matplotlib
matplotlib.use("qt5agg")
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
# from cartogrampy.components.cartogram import Cartogram

ldn_boro = gpd.read_file('data/LDN_Boro.geojson').to_crs(epsg=27700)
ldn_pop  = pd.read_excel('data/Pandas_Lon_Pop.xlsx')
data = ldn_boro.merge(ldn_pop, how='left', left_on='GSS_CODE', right_on='New Code')
c = Cartogram(data, value_field=2015 )
dorling = c.dorling(iterations=200,verbose=True)
ncc = c.noncont()

# Map the outcome of dorling
fig, ax = plt.subplots(1,2,figsize=(16,8))
data.plot(column = 2015, cmap = 'PuRd', scheme = 'Quantiles', k=4,legend=True, ax=ax[0])
dorling.plot(color='r',ax=ax[1])
ax[0].axis('equal')
ax[1].axis('equal')
plt.show()

# Map the outcome of NCC
fig, ax = plt.subplots(1,figsize=(8,8))
data.plot(color='black',ax=ax,alpha = 0.8,zorder=0)
ncc.plot(color='r',ax=ax)
ax.axis('equal')
ax.axis('equal')
plt.show()
