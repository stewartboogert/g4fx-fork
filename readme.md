# g4fx


## Getting Started 

### Build Instructions


Following assumes ubuntu based distros

Required packages

```
apt-get install docker \ 
                build-essential \
                libvtk9-qt-dev \
                libvtk9.1-qt \
                libcgal-dev \
                libgmp-dev \
                libmpfr-dev \
                libboost-all-dev \
                libxerces-c-dev \
                cmake \
                cmake-curses-gui \
                libxmu-dev \
                x11vnc \
                wget \
                git \
                libxi-dev \
                libglfw3-dev
```

For OpenCascade Support
```
apt-get install libocct-data-exchange-dev \
                libocct-draw-dev \
                libocct-foundation-dev \
                libocct-modeling-algorithms-dev \
                libocct-modeling-data-dev \
                libocct-ocaf-dev \
                libocct-visualization-dev \
                occt-draw \
                occt-misc
```


```
export G3FX_ROOT=git clone location
```

* [Building the Docker Container](docker/readme.md)
* [Building 3rd party libraries](3rdparty/readme.md) 

#### Build clhep

```
wget https://proj-clhep.web.cern.ch/proj-clhep/dist1/clhep-2.4.6.0.tgz && \
    tar -xf clhep-2.4.6.0.tgz && rm -rf clhep-2.4.6.0.tgz && \
    mkdir clhep-build && cd clhep-build && cmake -DCMAKE_INSTALL_PREFIX=$G4FX_ROOT/.install ../2.4.6.0/CLHEP/ && \
    make -j4 && make install && cd ../ && rm -rfv 2.4.6.0 clhep-build

```

#### Build Geant4
```
git clone https://gitlab.cern.ch/geant4/geant4.git &&\
    cd geant4 && git checkout geant4-11.1-release && cd ../ &&\
    mkdir geant4-build && cd geant4-build &&\
    cmake ../geant4 -DCMAKE_INSTALL_PREFIX=$G4FX_ROOT/.install -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_GDML=ON -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_USE_RAYTRACER_X11=ON -DGEANT4_USE_VTK=ON && make -j4 && make install && cd ../ &&\
    echo 'source /usr/local/bin/geant4.sh' >> ~/.bashrc

```

#### Building g4fx

```
mkdir $G4FX_ROOT/.build
cd $G4FX_ROOT/.build
cmake .. -DCMAKE_PREFIX_PATH=$G4FX_ROOT/.install -DCMAKE_MODULE_PATH=$G4FX_ROOT/.install/lib/cmake/Geant4
make -j8
```

#### Setting up Clion

Settings -> Build and Deploy -> CMake -> CMake Options

```
-DCMAKE_MODULE_PATH=$G4FX_ROOT/.install/lib/cmake/Geant4 -DCMAKE_PREFIX_PATH=$G4FX_ROOT/.install
```