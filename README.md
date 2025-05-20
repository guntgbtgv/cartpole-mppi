# Model Predictive Path Integral (MPPI) controller for cartpole
minimal c++ implementation of MPPI controller with mujoco 2.3.7

[![Watch the video](https://img.youtube.com/vi/-gL8eOqGG-k/default.jpg)](https://www.youtube.com/watch?v=-gL8eOqGG-k)

## How to use

`mujoco 2.3.7` need to be copied into the project 

```shell
git clone git@github.com:guntgbtgv/cartpole-mppi.git
cd cartpole-mppi
mkdir mujoco
cp <mujoco 2.3.7 dir> <cartpole-mppi/mujoco>
```

build the project and run

```shell
mkdir build
cd build
cmake ..
make
./mppi_test
```

## Reference

- [dm-control](https://github.com/google-deepmind/dm_control/tree/main)
- Williams, Grady, Andrew Aldrich, and Evangelos A. Theodorou. "Model predictive path integral control: From theory to parallel computation." Journal of Guidance, Control, and Dynamics 40.2 (2017): 344-357.