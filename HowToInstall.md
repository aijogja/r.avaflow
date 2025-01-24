
Install the grass.install.sh script manually, by open the file and install one by one

install the r.avaflow module

```
g.extension extension=r.avaflow.40G url=/home/aijogja/devel/avaflow/r.avaflow.40G/
```

run the Frankslide

```
r.in.gdal -o input=TIFF/fs_elev.tif output=fs_elev
g.region -p
g.region -s rast=fs_elev
g.region -p

```

```
installed.packages()[,1]
install.packages("terra")
install.packages("terra", dependencies=TRUE)
```
