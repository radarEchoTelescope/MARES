# **MARES Installation Instructions**


This document is intended to serve as a guide to installation of the MARES (a Macroscopic Approach to the Radar Echo Scatter) C++ code package, developed by Enrique Huesca Santiago. MARES is a semi-analytic, macroscopic code that calculates the radar echo signal from high-energy particle cascades in dense media. The motivation and physics behind the code can be found in \[1\], and in further detail in \[2\]. 

The first section details the installation of MARES on a computing cluster, and is where the full instructions can be found. MARES can also be run on a local system, and later sections will go over the necessary steps for this. The latter is slightly more complex and may require some fiddling to get working properly \- therefore, it is assumed that readers of this section will have looked over the first section and have some understanding of the code and its dependencies. 

This guide assumes some familiarity with SSH and linux commands \- [see here](https://www.star.bris.ac.uk/linuxtut/index.html) for a reference guide to get started with using these commands. It also currently requires access to the [MARES GitHub repository](https://github.com/radarEchoTelescope/MARES). 

# Contents

[MARES installation on a cluster](#mares-installation-on-a-cluster)

[Quick installation](#quick-installation)

[Full Instructions](#full-instructions)

[Usage](#usage)

[MARES installation on local systems](#mares-installation-on-local-systems)

[References](#references)

[Contact](#contact)

# MARES installation on a cluster

In this section, the instructions for installing MARES on a cluster are detailed. This is the easiest and quickest way to run MARES. In the beginning of this section, only the necessary terminal commands are given; primarily for those who have downloaded MARES multiple times but keep forgetting the minutiae of the install process. First-time users are recommended to look at the second part of this section, where the instructions have been written out in full in order to make clear what each step is doing and how it should be run. Following this, a very brief usage guide is provided, to test the installation. 

Terminal commands are signified with e.g. ```./configure``` and commands that should be adjusted for the user's choice of directory for the MARES install are signified with brackets as follows `<...>`.   

Before running any of the commands below, it is necessary to choose the correct CMakeLists.txt file. 
1. External Cluster: rename the CMake file `CMakeLists_cluster.txt` to `CMakeLists.txt`
2. Local MacOs system: rename the CMake file `CMakeLists_MacOS.txt` to `CMakeLists.txt`
3. Local Linux system: rename the CMake file `CMakeLists_Linux.txt` to `CMakeLists.txt`

## Quick installation

Here, only the terminal commands needed for setup are given, for quick reference. It is assumed that the user has set up the [necessary SSH key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent) for the github commands to work successfully. If not, replace git clone lines with an alternative method of sourcing the necessary packages from the repositories. 

Lines that should be adjusted for the user's choice of directories and path setup are signified with ```<...>```. When cloned from the repository, a new directory called MARES will be made, containing the code package. This code will create the MARES directory within a new directory of choice (MARESfolder) \- if this is not wanted, simply skip the first two lines. 


Make the folder containing the MARES software:

```
mkdir <... /MARESfolder>   
cd <... /MARESfolder>  
```
Clone the project:
```
git clone git@github.com:radarEchoTelescope/MARES.git 
````

Clone libconfig from the GitHub repository in the correct folder in the MARES project:
```
mkdir MARES/extern  
cd MARES/extern 
git clone git@github.com:hyperrealm/libconfig.git
```
Configure and install libconfig:
```
mkdir -p libconfig/bin
cd extern/libconfig   
autoreconf configure.ac
```
If the following error occures after the last terminal command: 
````
error: possibly undefined macro: AC_CHECK_INCLUDES_DEFAULT If this token and others are legitimate, please use m4_pattern_allow.
````
replace the `AC_CHECK_INCLUDES_DEFAULT`macro by `AC_HEADER_STDC` in the `configure.ac` file.

Next, configure libconfig:

```
./configure --prefix=<.../MARESfolder/MARES/extern/libconfig/bin>  
make  
make check  
make install   
make clean  
```
Compile project:
```
./INSTALL.sh 1
```

A successful completion of the INSTALL.sh script will signify that the code is successfully installed and ready to use.

A more detailed step-by-step guide can be found below.


## Full Instructions

1. ***Source the MARES package:***  
   First, choose the location of your new MARES install (e.g. `/user/software/`) and cd to that location. Then go to the Github repository ([https://github.com/radarEchoTelescope/MARES](https://github.com/radarEchoTelescope/MARES)) and retrieve the code package. This can be done via SSH (which requires the setup of [Github SSH keys](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)) using the terminal command:  
   

    ```
    git clone git@github.com:radarEchoTelescope/MARES.git
    ```

    Alternatively, the package can be [downloaded as a zip file directly from the repository](https://github.com/radarEchoTelescope/MARES/archive/refs/heads/main.zip). The zip file should then be uploaded to the cluster in the chosen install location, and expanded.  
     
     
1. ***Source and configure the libconfig package:***

   ***3a.*** We require the libconfig package to allow for the usage of config files. The full instructions can be found in the INSTALL file in the libconfig github ([https://github.com/hyperrealm/libconfig](https://github.com/hyperrealm/libconfig)), from which we have pulled the necessary sections here. 

   In this case, we first have to clone the libconfig repository from the GitHub page in the /`MARES/extern/` directory using the command
   
   ```
   git clone git@github.com:hyperrealm/libconfig.git
   ```
   Which will produce a directory, `libconfig`. 

   ***3b.***  Cd to the new `libconfig` directory, and make a new `bin` directory. The path should follow `/MARES/extern/libconfig/bin`.

   

   ***3c.*** We now have to configure and install libconfig. Cd back to `/MARES/extern/libconfig`, and run the command to produce the configure file:
   ```
   autoreconf configure.ac
   ```
   If the following error occures after the last terminal command: 
   ```
   error: possibly undefined macro: AC_CHECK_INCLUDES_DEFAULT If this token and others are legitimate, please use m4_pattern_allow.
   ```
   replace the `AC_CHECK_INCLUDES_DEFAULT`macro by `AC_HEADER_STDC` in the `configure.ac` file. This can occur if the version of *autoconf* is 2.69 or lower.

   Next, we configure the project:

    ```
    ./configure prefix=<.../MARES/extern/libconfig/bin>
    ```

   

   The `<.../MARES/extern/libconfig/bin>` parameter should be adjusted to follow the path of the user's MARES installation. Next, run a series of commands which compile, check and install the libconfig package respectively:
    ```
    make
    make check
    make install 
    ```
   If all output looks mostly healthy (there may be some minor errors, which can generally be ignored as any fatal errors will reveal themselves later) finish the installation with: 
    ```
    make clean 
    ```
   This removes unnecessary files from the source code directory.  

   Within the `/MARES/extern/libconfig/bin` directory, there should now be three new directories (include, lib, share).
   

2. ***Compile MARES***  
   Now all necessary packages have been added to make MARES run, we should be able to compile the code. Return to the `/MARES/` directory, where an `INSTALL.sh` script can be found.   
   
   To compile the project for the first time, run the script as follows:   
     
   	```
    ./INSTALL.sh 1
    ```  
     
   The '1' flag after the execution command tells MARES that it should start the compilation from scratch \- it will (re)make the necessary build directory and execution scripts for the code, in this case for the first time. It will also pull the IceRayTracing module from the corresponding GitHub page. This also clears any existing build and bin directories, completely starting over, so should mainly be used in the case something goes horribly wrong; any files in an existing build directory will be deleted.   
     
   The script will produce a lot of terminal output that should end with the lines:  
    ```  
    [95%] Building CXX object CMakeFiles/MARES.dir/MARES.cpp.o  
    [100%] Linking CXX executable MARES  
    [100%] Built target MARES  
    ```
   If this is the case, MARES can now be used.   
     
   

Note \- it is important that the paths for all the separate packages/installs are correct, as otherwise some very confusing libconfig errors can be produced, which are difficult to correct and make sense of. If any CMake errors occur after the ```./INSTALL.sh``` command, **the most likely error is that the paths are incorrect**.

# 

## Usage

Successful compilation of MARES will lead to a populated build and bin directory, where the executables for the different MARES simulation options can be found. 

Running and steering a MARES simulation requires a config file, of which an example has been provided in `/examples/example.cfg`. Here, the parameters of the simulation can be set by the user and provided to the code. Additionally, the type of output produced and the chosen output directory can be specified.  

When the config file is set up as desired, Move to the `MARES`directory and run: 
```
./bin/MARES examples/example.cfg
```
where the output of the simulation will be saved in the `output`directory. There should be some terminal output declaring that the config file has been parsed successfully (or not…), after which the code should execute and produce the chosen output.

To change parameters of the simulation, edit the contents of the example.cfg file. Lots of different variables to control the simulation, as well as options to change the ID of the file output and other parameters can be set by the user. 

Note \- if the source .cc, .cpp or .hh MARES scripts are modified, the code should be compiled every time before usage (or the changes will not be applied). In this scenario, the command is: 
```
./INSTALL.sh 0
```
When run with the `0` flag, the `INSTALL.sh` script will only compile the scripts, and any other files that were changed.

In order to update the external libraries (libconfig and IceRayTracing) from their github page, the command is: 
```
./INSTALL.sh 2
```

# MARES installation on local systems

Here, we detail the necessary extra steps needed to get MARES working on local systems. This was primarily worked out for Apple MacOS (ARM64), which can be a little fiddly, but has also been shown to work on a Linux system. The method has only been tested twice (for these respective systems), and it is possible more modifications will be required for different systems. A local windows installation is not included in this manual. We recommend installing Windows Subsystem for Linux (WSL) and following the installation guide for MARES on *Linux* outlined below. WSL can be installed following the [guide](https://learn.microsoft.com/en-us/windows/wsl/install). 

First, step 1 from the full installation instructions can be followed, providing us with the MARES in the correct location. 

Following this, it is likely that a few extra packages will be required: CMake, GSL, PkgConfig and LibConfig. Depending on the local system and user, some of these may already be installed, in which case they can be skipped (they are generally installed globally on computing clusters, which is why this step has been skipped in earlier sections. 

1. ***MacOS installation***  
   These packages can be retrieved (and configured) via a package manager \- for MacOS, we recommend [homebrew](https://brew.sh/).   
     
   Using the package manager of choice (here, homebrew), go ahead and install the necessary packages:   
   [PkgConfig](https://formulae.brew.sh/formula/pkgconf#default)  
   [GSL](https://formulae.brew.sh/formula/gsl#default)  
   [CMake](https://formulae.brew.sh/formula/cmake#default)   
   [LibConfig](https://formulae.brew.sh/formula/libconfig#default)  
     
   With these packages installed, make sure to have renamed the `CMakeLists_MacOS.txt` to `CMakeLists.txt`. This is important for the local system paths/packages, especially for libconfig. The CMakeLists.txt file included, was tested for a M4 MacOS system (Sequoia 15.5) and M5 MacOS system (Tahoe 26.2).
     
   If the MARES compilation (`./INSTALL 1`) finishes without errors, then the code should be successfully installed on the local system.   
     
2. ***Linux installation***  
   These packages can be retrieved (and configured) via a package manager \- for Linux, we recommend [apt](https://wiki.debian.org/Apt). Note that the manager [snap](https://snapcraft.io/) can be used also, but it has not been tested by the authors. Potential errors can arise due to its use of application containers.  
     
   First, update the linux system by running the following two commands:  
   ```   
   sudo apt update  
   sudo apt upgrade
   ```
     
   Next, install **GSL** and **CMake**. It is noted that the GSL package is included in the gcc compiler for C++.  Hence, the command below installs the compiler \- which is also needed for running MARES.
   ```
   sudo apt install build-essentials
   sudo apt install cmake
   ```
   

   Finally, **PkgConfig** and **Libconfig** are installed. The former can be installed using
   ```
   sudo apt install pkg-config
   ```
   

   For unknown reasons, a simple ```sudo apt install libconfig …``` command retrieves an old version of libconfig (1.5) that is not compatible with MARES. Instead, the necessary debian packages need to be downloaded from the repository and installed manually.   
   Retrieve  *libconfig++-dev\_1.7.3-2\_amd64* from [here](https://launchpad.net/ubuntu/questing/+package/libconfig++-dev) and *libconfig++11\_1.7.3-2\_amd64* from [here](https://launchpad.net/ubuntu/questing/+package/libconfig++11). To install the packages, run the commands below with the correct path to your download folder. 
   ```
   sudo apt install  <.../libconfig++-dev_1.7.3-2_amd64>  

   sudo apt install  <.../libconfig++11_1.7.3-2_amd64>
   ```

   

	  
    With these packages installed, the only task remaining is that the CMakeLists.txt file (most likely) needs to be modified in order to account for the local system paths/packages, especially for libconfig. This is the tricky part, as the required modifications seem to differ depending on the system. An example working CMakeLists.txt file has been included below, for a Ubuntu 24.04 system, which may be a good place to start from. 

    With these packages installed, make sure to have renamed the `CMakeLists_Linux.txt` to `CMakeLists.txt`. This is important for the local system paths/packages, especially for libconfig. The CMakeLists.txt file included, was tested for a Ubuntu 24.04 system.

    If the MARES compilation (`./INSTALL 1`) finishes without errors, then the code should be successfully installed on the local system.   

3. ***Windows installation***  
   See *Linux installation* using the WSL software on Windows. 

	


# Footnotes:

[^1]:  If the MARES package retrieved from the github already has existing directories for some external packages (IceRayTracing or libconfig) it is recommended to delete these and start with fresh installs.  

[^2]:  Alternatively, if the `…/MARES/extern` folder does not exist in the MARES install, first make a new `extern` folder in the MARES directory, then follow the subsequent instructions.

[^3]:  Technically, we only require this version of libconfig because the configuration instructions we are following have yet to be updated to reflect more recent versions of the code. MARES runs perfectly well for newer libconfig versions, as shown for local installations, so this section may be modified to reflect that, as soon as the authors can pin down how the newer libconfig should be properly configured on a cluster. 

 # References

\[1\] E. Huesca Santiago et al., (The Radar Echo Telescope Collaboration), '*Macroscopic approach to the radar echo scatter from high-energy particle cascades*', [Phys. Rev. D. **109** 083012 (2024)](http://doi/10.1103/PhysRevD.109.083012)  
\[2\] E. Huesca Santiago, '*Understanding Radar Echoes from High-Energy Particle Cascades*', PhD Thesis, Vrije Universiteit Brussel (2024)

# Contact
**Krijn de Vries**: krijndevries AT gmail.com

**Enrique Huesca Santiago**: ehuescasantiago AT gmail.com