# Head-Neck Cancer Project

## Background

The head and neck region, including tumorous volumes within it, is not anatomically rigid, but is very radiosensitive. This means that its 3D shape is in constant flux during radiotherapy, and is impacted by tumour shrinkage, inflammation, weight loss, and changes that affect the patient’s muscle density, fat distribution, edema, and fluid accumulation. Indeed, HNC patients who undergo radiotherapy often experience treatment-related toxicities that result in anatomical changes. These anatomical changes in turn can render the treatment plan invalid, requiring replanning, which if performed ad hoc, as is currently the case, can be disruptive for the treatment planning team. Therefore, in this research project we are investigating if we can predict if and when a head and neck cancer patient requires replanning using an AI algorithm.

## Goals

The ultimate goal of this project is to predict when a HNC patient would require replanning or if they would not require radiotherapy replanning, at the early stages of the treatment so that there are no last minute changes and the treament follows a planned procedure.



## Folders

* [Spactial-feature extractor](/HNC_project/head-and-neck_spatial-feature_extractor-main/) : Code to extract spatial features from CT and CBCT images of head and neck cancer patients.

* [ Neck-area calculation](/HNC_project/Neck_area_calculation/) : Code was written to calculate the area/volume of a region of the neck using DICOM files from a Radiotherapy treatment

* [PTV_body_distance-main](/HNC_project/PTV_body_distance-main/) : This folder contains scripts to calculate and analyse the minimum distance between the PTV and the body.

* [ Data Analysis and Machine Learning ](/HNC_project/Data_Analysis_and_Machine_Learning/) : This folder contains python code used to combine, analyse and visualise HNC patient data. It also contains machine learning code used to predict whether patients would require radiotherapy replanning or not.
