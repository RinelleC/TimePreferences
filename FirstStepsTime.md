# First Steps for Time Preferences Analysis

Rinelle Chetty
28 August 2025

## Create Time Preferences folder
mkdir TimePreferences 
cd TimePreferences
touch FirstStepsTime.md      <-- the current file 

## Create folders to use during analyses in Stata 
cd TimePreferences 
mkdir dofiles
mkdir logfiles 
mkdir estimates 
mkdir stata_tables 
mkdir stata_graphs 
mkdir figures 


## Add submodules 

### First submodule 
git clone https://github.com/andrehofmeyr/COVID-19-data.git
cd COVID-19-data 
git submodule init 
git submodule update --recursive 

### Second submodule 
../ 
git clone https://github.com/andrehofmeyr/covid_clean_rsa.git
cd covid_clean_rsa 
git submodule init 
git submodule update --recursive 

### Third submodule (the parent folder)
../
git clone https://github.com/andrehofmeyr/covid_data_rsa.git 
cd covid_data_rsa 
git submodule init 
git submodule update --recursive 


## Run Clean.do from covid_clean_rsa
Note that this has been incorporated into my first "cleaning" do file. 

## Run Main.do from covid_data_rsa
My analyses have been conducted across multiple do files. 


## Add this local folder to Github
### In VS Code: 
git init 
git add . 
git commit -m "first commit" 

### In Github online: 
1. Login.
2. Create new repository with the same name as the local folder. 
3. Copy the link created. 

### Return to VS Code and in the terminal type: 
git remote add origin <<link>> 
(i.e., git remote add origin https://github.com/RinelleC/TimePreferences.git )

git push -u origin master

### Refresh the folder in GitHub online 

## Now run the do files

END.