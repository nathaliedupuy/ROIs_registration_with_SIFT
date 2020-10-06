###########################################################
# INSTALL ROIS REGISTRATION WITH SIFT (OPENCV) ON WINDOWS #
###########################################################
CONTENT
--------------------------------------------------------------------------------
	INSTALLATION STEPS - PART I  - PREP THE VIRTUAL ENV
	INSTALLATION STEPS - PART II - OPENCV WITH SIFT
	RUN PYTHON IN MATLAB WITH ANACONDA (WINDOWS)
	TROUBLESHOOT


# This document is long because it is step-by-step, but the installation is not that complicated...

# INSTALLATION STEPS - PART I - PREP THE VIRTUAL ENV
# ------------------------------------------------------------------------------


**Important before you start - if you do not use Anaconda**

It is recommended to use Anaconda if you are not familiar with python package management.
Most of the steps are descibed for anaconda, but anaconda is not required and
you can still follow the guide using pip.


**Important before you start - if you use Anaconda**

If you have both Anaconda2 and Anaconda3 this will likely cause problems.
It is best to uninstall everything and install only Anaconda3. See Troubleshoot.


**(WINDOWS) Open Anaconda command prompt / or Windows cmd if you do not use Anaconda**

When you install Anaconda you may remember that the box to *add Anaconda to the PATH*
was unchecked by default, which means that conda commands will not be recognized
if you try to use them in the default Windows command prompt.
=> It is recommended to keep it that way (*add Anaconda to the PATH* unchecked).

So instead we use Anaconda prompt.

Some basic commands that may help (work for both):

	- To navigate, use :
		cd thepath
	- To back one step:
		cd ..
	- To change disk, e.g. from C:\ to D:\, simply write
		d:
	- To display directory content:
		dir
	- To start/exit python:
    	python
    	exit() or quit()


**Create a virtual environment**

If you do not wish to use anaconda, use virtualenv
    https://virtualenv.pypa.io/en/latest/
    and then Pip to manage the packages


With Anaconda

	conda create --name name_of_venv python=x.xx

A new directory named *name_of_venv* will be created by default in the
/envs/ where Anaconda3 is (e.g. C:\Users\the_user\Anaconda3\envs). 
Choose a useful name for your virtual environment.
This environment uses a specific version of Python x.xx.
If you do not know which one to pick, set to python3.6.

Some basic commands for the virtual env that may help:

	- List of venvs: conda env list
	- Activate environment: conda activate *name_of_venv*
	- Deactivate environment: conda deactivate
	- Delete environment:  conda remove --name *name_of_venv* --all	



**Packages installation in the virtual environment**

Activate environment:   

	conda activate name_of_venv

If successful, you will see the start of your command line changed with:
    (*name_of_venv*) path\where\you\are>

In the virtual environment, install the following packages:
	
	conda install numpy
	conda install pyyaml
	conda install scipy

Some basic commands for the virtual env that may help:

	- Listing all packages in the current virtual environment when the venv is activated:
    	conda list
    - To check a specific package in the current virtual environment:
    	conda list package_name
    - To remove a specific package in the current virtual environment:
    	conda remove package_name


**Check everything is ok so far**

Start MATLAB from the virtual env (simply type 'matlab'). Navigate to the 'test_install'
folder in the regRois folder (no need to add anything to the path). 
You should find two dummy python files there.
In MATLAB, check that you have no error when you run:

	system('python dummy1.py', '-echo');

If you have a DLL error, uninstall scipy and re-install using pip:
	pip install scipy



# INSTALLATION STEPS - PART II - OPENCV WITH SIFT
# ------------------------------------------------------------------------------

We need to build OpenCV from source because of the SIFT library which is now an 
extra feature, as it is no longer open access except for research.


**Visual Studio**

If you do not have it Download & install
	https://visualstudio.microsoft.com/
	Choose free version Visual Studio Community.
	
When you install Visual Studio, make sure to include the module called 
	"Desktop development with C++"  (it is needed for building)
	Make sure to add the SDK, and restart your computer if you've made any changes.
	If you do not see it, make sure you have a recent version of VS.
If you missed it, no worries you can install it later:
	simply launch Visual Studio Installer, then click modify and select the module.
	

**Download CMake** 

https://cmake.org/download/
Look for Binary distributions - and choose Windows winXX-xXX *ZIP*


**Download OpenCV from Github**
		
https://github.com/opencv/opencv
Also Download & Extract OpenCV-contrib from Github (this is where SIFT is)
https://github.com/opencv/opencv_contrib

In the OpenCV-contrib/modules folder you can delete all folders
but 'xfeatures2d' which is the one with SIFT (unless you need extra
modules of course!). This will make the installation lighter and 
avoid unnecessary potential issues...

	- In the OpenCV folder, create a new folder called 'build'.

*Important*
If at any point something fails during the steps below, 
delete the content of the 'build' folder before re-do the next steps.



**CMake in action**
	
	- Back in the command prompt, deactivated the virtual env (conda deactivate)
		Check you are back to base.

It would have been great to use the gui, but things are never that
straightforward.
	
	- In the command prompt, navigate to the folder which contains cmake.exe
		The exe files should be in the_cmake_dir/bin/
		
	- Replace the parts in brackets [] (and remove the brackets!) 
		in the lines below with your paths
		(attention for DPYTHON3_LIBRARY you also need to specify the python version) 
		and then copy the block to command line.

*VERY VERY Important*
1) Do not use backslash '\' but slash '/' for the paths.
2) NO SPACE after the  '=', directly insert the paths/values. 
Same goes for -S and -B. Directly insert the paths of source and build.
(Note that the '^' indicate line continuation)

	cmake -DCMAKE_BUILD_TYPE=RELEASE ^
		-DINSTALL_PYTHON_EXAMPLES=OFF ^
		-DINSTALL_C_EXAMPLES=OFF ^
		-DBUILD_EXAMPLES=OFF ^
		-DBUILD_opencv_python3=ON ^
		-DBUILD_opencv_python2=OFF ^
		-DPYTHON_DEFAULT_EXECUTABLE=[ path/to/your/virtualenv ]/python.exe ^
		-DPYTHON3_EXECUTABLE=[ path/to/your/virtualenv ]/python.exe ^
		-DPYTHON3_INCLUDE_DIR=[ path/to/your/virtualenv ]/include ^
		-DPYTHON3_NUMPY_INCLUDE_DIRS=[ path/to/your/virtualenv ]/Lib/site-packages/numpy/core/include ^
		-DPYTHON3_LIBRARY=[ path/to/your/virtualenv ]/libs/python[ XX ].lib ^
		-DPYTHON3_PACKAGES_PATH=[ path/to/your/virtualenv ]/Lib/site-packages ^
		-DOPENCV_ENABLE_NONFREE=ON ^
		-DOPENCV_EXTRA_MODULES_PATH=[ path/to/OpenCV-contrib ]/modules ^
		-S[ path/to/OpenCV ] -B[ path/to/OpenCV ]/build

Note that the line continuation ^ might not work... 
in that case remove them and write everything in one single line.

Check the output before continuing:

	- Check the Python3 configurations are correctly set with your virtual environment. 
		You should see:
			Interpreter, Libraries, numpy and install path.
	- In the section
			OpenCV modules:
				To be built:
		you should see at the end 'xfeatures2d' at the end of the list.


**Build with Visual Studio**

In your explorer, go to the OpenCV/build folder.
	
	- Open OpenCV.sln file with Visual Studio.
	- Check 'build mode' is 'Release' instead of 'Debug' (Top toolbar)
	- In the solution explorer, go to 'CMakeTargets' 
		right-click on ALL_BUILD and build it (this part takes some time)
	- Still in 'CMakeTargets' in the solution explorer, 
		right click on INSTALL and build it.


**Check installation of opencv**

*To check everything ok when calling from MATLAB:*

In the anaconda prompt, activate the virtual environment and start MATLAB (simply type 'matlab'). 
Navigate to the 'test_install' folder in the regRois folder (no need to add anything to the path). 
You should find two dummy python files there.
In MATLAB, check that you have no error when you run:

	system('python dummy2.py', '-echo');

*If you wish to test with python only and not with MATLAB:*

	- Activate the virtual environment.
	- Start python (simply type 'python').
	- Then check opencv is installed correctly:
		import cv2
		cv2.__version__
	- Initiate a dummy SIFT detector to check the extra module was installed:
		sift = cv2.xfeatures2d.SIFT_create()

Note: if you start python from base, you should not be able to import
opencv (No module named 'cv2') as we installed it in the virtual env.



**Free up space!**

If you need to free up space, you can safely delete the following files:

	- opencv_contrib folder
	- in opencv folder, keep the BUILD folder, LICENCE, README.md and SECURITY.md
		and delete the rest
	- in the build folder only keep the install folder; delete the rest


<<< More info >>>
- opencv install on windows
https://docs.opencv.org/master/d5/de5/tutorial_py_setup_in_windows.html



# RUN PYTHON IN MATLAB WITH ANACONDA (WINDOWS)
# ------------------------------------------------------------------------------

If you need to call python scripts that require specific packages installed with 
Anaconda you have to start MATLAB from the Anaconda Prompt.

In the Command Line simply type 'matlab'.

*Virtual Environments*

If you want to run Python from a virtual environment you need to specify
the full path to the python binary. For example:

	system('/path/to/venv/python myscript.py')


Alternatively, you can activate the virtual environment in the cmmand prompt
and call matlab from within. This is what we did to test everything was working 
during the installation. The libraries paths will be automatically set to the 
virtual environment (you can check by typing 'pyversion' in MATLAB).


# TROUBLESHOOT
# ------------------------------------------------------------------------------

**Anaconda command did not work**

Double check the version of anaconda you are using.
Some methods differ slightly (e.g. 'conda activate' vs 'activate' alone...)


**Anaconda2 + Anaconda3**

This likely caused problems for virtual env etc...
Make sure everything important is saved somewhere (e.g. any created virtual env).
Uninstall all, and reinstall the latest (Anaconda3).
Note that it will install python3, but you can always create a virtual
environment with python2 to work with older projects.

For more info, see:
https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#viewing-a-list-of-your-environments


**Warning during cmake**

Just ignore them...

**Error during opencv prep installation with cmake**
	
	- Check paths: on WINDOWS paths are with backslash '\' or slash '/'
		However for CMake use only slash '/'
		(note on Linux use slash '/')
	- Make sure you downloaded the same version for opencv and opencv_contrib


**Visual Studio asks for a licence**

It happens... It is free but you need to sign-in...


**Error DLL when running python from MATLAB**

Make sure you are calling the python from the right virtual env.

(WINDOWS) If you use Anaconda, first make sure you are running matlab from the anaconda prompt, 
and that you call python correctly (see section RUN PYTHON IN MATLAB WITH ANACONDA (WINDOWS)).

*Issues with scipy* 
I got issues when importing scipy.stats in matlab, 
while it worked fine when calling python from cmd line.
A workaround is to uninstall scipy and re-install using pip:
	pip install scipy


*For other packages* 
Do not delete anything in the opencv folder (skip the 'free-some-space' step)
and try again. If it fails again, re-do the installation...


