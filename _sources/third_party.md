# List of provided meshes from third party

## LifeX

The LifeX example meshes can be found at https://zenodo.org/record/5810269#.YeEjWi8w1B0, which also contains a DOI: https://doi.org/10.5281/zenodo.5810269.


## Computational meshes of the cardiac left ventricle of 50 heart failure subjects

Collection of 50 anatomies of the left ventricle of a cohort of subjects selected for cardiac resynchronisation therapy, togheter with their response data.Meshes: provided in VTK format for easier compatibility (original work done with cubic Hermite EX format meshes)Response: excel spreadsheet with the case ID and the class (responder or not responder). "data process.gif" provides an overview of the segmentation and mesh personalization process. Panel A: semi-automatic segmentation by 3D extrapolation (yellow surface and contours) of 2D segmentation contours (red contours and projections) using an open source platform (29). Panel B: mesh template. Panel C: resulting mesh (white surface) overlaid with the segmentation surface colour coded by the distance between them (jet colour map, from 0mm in blue to 1mm in red).

https://figshare.com/articles/dataset/Computational_meshes_of_the_cardiac_left_ventricle_of_50_heart_failure_subjects/5853948/1

## Virtual cohort of adult healthy four-chamber heart meshes from CT images

Rodero, Cristobal; Strocchi, Marina; Marciniak, Maciej; Longobardi, Stefano; Whitaker, John; O'Neill, Mark D.; Gillette, Karli; Augustin, Christoph; Plank, Gernot; Vigmond, Edward J.; Lamata, Pablo; Niederer, Steven A.

Dataset Description: We present the first database of four-chamber healthy heart models suitable for electro-mechanical (EM) simulations. Our database consists of twenty four-chamber heart models generated from end-diastolic CT acquired from patients who went to the emergency room with acute chest pains. Since no cardiac conditions were detected in follow-up, these patients were taken as representative of "healthy" (or asymptomatic) hearts. These meshes were used for EM simulations and to build a statistical shape model (SSM). The output of the simulations and the weights of the SSM are also provided.

Cardiac meshes: We segmented end-diastolic CT. The segmentation was then upsampled and smoothed. The final multi-label segmentation was used to generate a tetrahedral mesh. The resulting meshes had an average edge length of 1 mm. The elements of all the twenty meshes are labelled as follows:

https://zenodo.org/record/4590294#.YfbZIS8w1B0

## Virtual cohort of extreme and average four-chamber heart meshes from statistical shape model

Rodero, Cristobal; Strocchi, Marina; Marciniak, Maciej; Longobardi, Stefano; Whitaker, John; O'Neill, Mark D.; Gillette, Karli; Augustin, Christoph; Plank, Gernot; Vigmond, Edward J.; Lamata, Pablo; Niederer, Steven A.

Dataset Description: We present a database of four-chamber heart models derived from a statistical shape model (SSM) suitable for electro-mechanical (EM) simulations. Our database consists of 39 four-chamber heart models generated from end-diastolic CT-derived meshes (available in the repository called ("Virtual cohort of adult healthy four-chamber heart meshes from CT images"). These meshes were used for EM simulations. The output of the simulations and the weights of the SSM are also provided.

Cardiac meshes: To build the SSM, we rigidly aligned the CT cohort and extracted the surfaces, representing them as deRham currents. The registration between meshes and computation of the average shape was done using a Large Deformation Diffeomorphic Metric Mapping method. The deformation functions depend on a set of uniformly distributed control points in which the shapes are embedded, and on the deformation vectors attached to these points. It is in this spatial field of deformation vectors (one per each control point) where the Principal Component Analysis is applied. Case #20 of the CT cohort was not included. More information on the details can be found in Supplement 3 of the reference paper. We created two extra cohorts by modifying the weight of the modes explaining 90%of the variance in shape (corresponding to modes 1 to 9). We created these meshes with either ±2 or ±3 standard deviations (SD) of each mode added to the average mesh (extreme2 and extreme3 cohorts respectively). We also created two additional meshes with ±1 SD for mode 2 (extreme1 cohort). The elements of all the meshes are labelled as follows:

https://zenodo.org/record/4593739#.YfbZoC8w1B0


## A Virtual Cohort of Twenty-four Left-ventricular Models of Ischemic Cardiomyopathy Patients

Caroline Mendonca Costa; Aurel Neic; Eric Kerfoot; Karli Gillette; Bradley Porter; Benjamin Sieniewicz; Justin Gould; Baldeep Sidhu; Zhong Chen; Mark Elliott; Vishal Mehta; Gernot Plank; Aldo Rinaldi; Martin Bishop; Steven Niederer

Description
Motivation: Computational models of the heart are increasingly being used in the development of devices, patient diagnosis and therapy guidance. While software techniques have been developed for simulating single hearts, there remain significant challenges in simulating cohorts of virtual hearts from multiple patients. Dataset Description: We present the first database of left-ventricular (LV) models suitable for electrophysiology simulations. Our database consists of twenty-four LV models including infarct scar morphology. These were generated from LGE-MRI acquired from ICM patients undergoing CRT. We used 24 image-based patient-specific models of LV anatomy and scar morphology. Briefly, LV endocardium and epicardium contours were manually drawn in each short-axis slice of LGE-MRI. Scar and BZ were segmented and reconstructed in 3D. A finite element tetrahedral mesh (mean edge length of 0.8mm) was generated and 3D reconstructed scar and BZ segmentations were mapped onto it. Rule-based fibres were assigned to the models. The elements of all the twenty-four meshes are labelled as follows: 1) Left ventricular myocardium 3) Scar core 4) Scar border zone We defined a system of universal ventricular coordinates on the meshes: an apico-basal coordinate varying continuously from 0 at the apex to 1 at the base; a transmural coordinate varying continuously from 0 at the endocardium to 1 at the epicardium; a rotational coordinate varying continuously from – π at the left ventricular free wall, 0 at the septum and then back to + π at the left ventricular free wall. We built the first cohort of twenty-four LV meshes from ICM patients LGE-MRI data. These geometries can be used for large cohort computational studies.

https://doi.org/10.18742/RDM01-570

## Data from: A Publicly Available Virtual Cohort of Four-chamber Heart Meshes for Cardiac Electro-mechanics Simulations

Strocchi, Marina; Augustin, Christoph M.; Gsell, Matthias A. F.; Karabelas, Elias; Neic, Aurel; Gillette, Karli; Razeghi, Orod; Prassl, Anton J.; Vigmond, Edward J.; Behar, Jonathan M.; Gould, Justin S.; Sidhu, Baldeep; Rinaldi, Christopher A.; Bishop, Martin J.; Plank, Gernot; Niederer, Steven A.

Description
Motivation: Computational models of the heart are increasingly being used in the development of devices, patient diagnosis and therapy guidance. While software techniques have been developed for simulating single hearts, there remain significant challenges in simulating cohorts of virtual hearts from multiple patients.

Dataset Description: We present the first database of four-chamber heart models suitable for electro-mechanical simulations. Our database consists of twenty-four four-chamber heart models generated from end-diastolic CT acquired from heart failure patients recruited for cardiac resynchronization therapy upgrade. We also provide a higher resolution version for each of the twenty-four meshes.


https://doi.org/10.5281/zenodo.3890034



## Virtual cohort of 1000 synthetic heart meshes from adult human healthy population

Rodero, Cristobal; Strocchi, Marina; Marciniak, Maciej; Longobardi, Stefano; Whitaker, John; O'Neill, Mark D.; Gillette, Karli; Augustin, Christoph; Plank, Gernot; Vigmond, Edward J.; Lamata, Pablo; NIederer, Steven

Dataset Description: We present a database of four-chamber heart models derived from a statistical shape model (SSM) suitable for electro-mechanical (EM) simulations. Our database consists of 1000 four-chamber heart models generated from end-diastolic CT-derived meshes (available in the repository called ("Virtual cohort of adult healthy four-chamber heart meshes from CT images"). These meshes were used for EM simulations. The weights of the SSM are also provided.

Cardiac meshes: To build the SSM, we rigidly aligned the CT cohort and extracted the surfaces, representing them asdeRham currents. The registration between meshes and computation of the average shape was done using a Large Deformation Diffeomorphic Metric Mapping method. The deformation functions depend on a set of uniformly distributed control points in which the shapes are embedded, and on the deformation vectors attached to these points. It is in this spatial field of deformation vectors (one per each control point) where the Principal Component Analysis (PCA) is applied. Case #20 of the CT cohort was not included. More information on the details can be found in Supplement 3 of the reference paper. We created this cohort by modifying the weight of the modes explaining 90%of the variance in shape (corresponding to modes 1 to 9) within 2 standard deviations (SD) of each mode added to the average mesh. The elements of all the meshes are labelled as follows:

https://doi.org/10.5281/zenodo.4506930


## Four-Chamber Human Heart Model for the Simulation of Cardiac Electrophysiology and Cardiac Mechanic

Tobias Gerach; Steffen Schuler; Andreas Wachter; Axel Loewe

This repository contains a four-chamber model of the human heart which is ready to use for simulations of cardiac electrophysiology and cardiac mechanics problems.

The cardiac anatomy was manually segmented from magnetic resonance imaging (MRI) data of a 33 year old male volunteer. The volunteer provided informed consent and the study was approved by the IRB of Heidelberg University Hospital (Fritz et al., 2014).
The MRI data were acquired using a 1.5 T MR tomography system and consist of a static whole heart image stack taken during diastasis as well as time-resolved images in several long and short axis slices. Based on the segmentation, we first labeled the atria and the ventricles. The geometry was extended by a representation of the mitral valve, the tricuspid valve, the aortic valve and the pulmonary valve. Additionally, we closed the endo- and epicardial surfaces of the atria and added truncated pulmonary veins, vena cavae as well as the ascending aorta and pulmonary artery. Furthermore, we added a concentric layer of tissue around the entire heart which phenomenologically represents the influence of the pericardium and the surrounding tissue.

Two tetrahedral meshes were created using Gmsh (Geuzaine et al., 2009): (1) the mechanical reference domain (M.vtu) with 128,976 elements (on average 3.17 mm edge length) and (2) the electrophysiological reference domain (EP.vtu) as a subset of M with 50,058,295 elements (on average 0.4 mm edge length).
We used rule-based methods to assign the myofiber orientation on EP: Wachter et al. (2015) was used for the atria and Bayer et al. (2012) for the ventricles. The fiber angle in the ventricles was chosen as +60° and -60° on the endocardial and epicardial surface, respectively. The sheet angle was set to -65° on the endocardium and 25° on the epicardium. Github repositories to these fiber generation tools are given in the sidebar. All geometry files are given in millimeter (mm).


https://doi.org/10.5281/zenodo.5573921


## Anatomically Detailed Human Atrial FE Meshes


Augustin, Chistoph M.; Fastl, Thomas E.; Neic, Aurel; Bellini, Chiara; Whitaker, John; Rajani, Ronak; O'Neill, Mark D.; Bishop, Martin J.; Plank, Gernot; Niederer, Steven A.

The left atrium (LA) has a complex anatomy with heterogeneous wall thickness and curvature. We include 3 patient-specific anatomical FE meshes with rule-based myofiber directions of each of the anatomies included in our study (The impact of wall thickness and curvature on wall stress in patient-specific electromechanical models of the left atrium, BMMB, 2020, https://pubmed.ncbi.nlm.nih.gov/31802292/).
Additionally we include
- a noised model with Gaussian noise added (mean 0 um , standard deviation 100 um ) to the initial geometry of patient case 3 and subsequently smoothed using ParaView; and
- a mesh with a constant wall thickness of 0.5 mm generated based on the endocardial surface of patient case 3.

The meshes are in the binary format for the Cardiac Arrhythmia Research Package simulator, see https://carpentry.medunigraz.at/carputils/index.html">https://carpentry.medunigraz.at/carputils/index.html and https://opencarp.org. For each of the geometries, we include a list of nodal coordinates (.bpts file), a list of triangular elements (.belem file), fiber fields (.blon file), surface files (*.surf files), and surface points (*.surf.vtx files).

Using the open source mesh manipulation utiliy "MeshTool" (https://bitbucket.org/aneic/meshtool/src/master/README.md">https://bitbucket.org/aneic/meshtool/src/master/README.md)
meshes can be converted to VTK or EnSight file formats.

https://doi.org/10.5281/zenodo.3843216


## Biventricular statistical shape model of the human heart adapted for computer simulations

Schuler, Steffen; Loewe, Axel


This is an adapted version of the biventricular statistical shape model from Bai et al. (2015) that can be used as a basis for computer simulations of cardiac electrophysiology or cardiac mechanics. The original model consists of disconnected surfaces of the left ventricular (LV) myocardium and the right ventricular (RV) blood pool. The surface of the RV blood pool was clipped at the base to add an orifice representing the RV in- and outlets. The resulting RV endocardial surface was shifted along its normals by a fixed wall thickness of 3 mm to obtain an RV epicardial surface. Then all surfaces were merged to form one closed surface of the biventricular myocardium. This surface was remeshed using Instant Meshes (Jakob et al., 2015) and tetrahedralized using Gmsh (Geuzaine. et al., 2009), resulting in the following two meshes of the mean shape:

https://zenodo.org/record/4419784#.Yfbc_S8w1B0


## Anatomy of the right ventricle of the heart in the young adult

Pablo Lamata; Adam J. Lewandowski; Afifah Mohamed

Anatomy of the right ventricle (RV) of 89 young subjects

Cohort 40 Preterm and 49 Term. Age: 18 - 40 years old (average Preterm: 22.7 years, Full term: 23.6 years) Data from YACHT: https://www.rdm.ox.ac.uk/about/our-clinical-facilities-and-mrc-units/cardiovascular-clinical-research-facility/ongoing-clinical-studies/yacht

Terms Please acknowledge this source of data, and associated publication:

A Mohamed et al. “Multimodality Imaging Demonstrates Reduced Right Ventricular Function Independent of Pulmonary Physiology in Moderately Preterm-Born Adults”

https://doi.org/10.6084/m9.figshare.11695236.v1


## Left ventricle of the maternal heart 9 years after pregnancy

Pablo Lamata


Description
Anatomy of the left ventricle LV of 153 mothers 5-10 years after pregnancy, with and without a history of hypertensionThe LV end diastolic anatomy was reconstructed from the contours using high order interpolation meshes as described in our previous work. The 153 meshes were then spatially aligned by their centre of mass, and by an orientation defined by their basal plane and the left to right direction set by the centre of mass of the LV and right ventricle. A statistical shape model was then built with a principal component analysis (PCA), finding the modes of anatomical variation of this cohort. Data shared containts:1. AtlasMeshData.zip, where each case has 1.1 binary.vtk: The mask of the left ventricular (LV) myocardium and right ventricular blood pool from the short axis magnetic resonance image in end diastole1.2 binary_mesh files (.exelem and .exnode): the cubic Hermite mesh that was personalised to the anatomy of the left ventricle1.3 Image_binary_mesh.jpg: image of the overlay between the fitted mesh and the mask of the LV myocardium1.4 Image_Initialization.jpg: image of the overlay between the template LV mesh used for personalization and the mask of the LV myocardium1.5 ImageDistances2ClosestPointbinary_meshxxxxxxxxxx.jpg: image of the distances between the LV mesh and the mask of the LV myocardium (xxxxxxxxxx is a digit encodign for the time when the image was generated)1.6 PersonalizationReportbinary.txt: file with metrics of mesh fitting accuracy and quality2. VTKMeshes.zip: the collection of all the LV shapes in VTK format, after correction of centre of mass and circumferential orientation. 3. AtlasModesPVSn43_rvspace_s4_c3.pdf: illustration of all the 10 first modes of anatomical variation, together with the box-plot of the differences between hypertensive pregnancy (blue) and normotensive pregnancy (green).4. Mode6revealsDifferences.png: screen capture of mode 6 from last file, revealing the largest differences between the groups of hypertensive and normotensive pregnancy

https://doi.org/10.6084/m9.figshare.11303306.v1


## Repository for modelling acute myocardial ischemia: simulation scripts and torso-heart mesh

Martinez-Navarro, Hector and Rodriguez, Blanca and Bueno-Orovio, Alfonso and Minchole, Ana. (2019).


Code and data required to run a biophysically detailed model of acute myocardial ischemia in humans on Chaste

https://ora.ox.ac.uk/objects/uuid:951b086c-c4ba-41ef-b967-c2106d87ee06
