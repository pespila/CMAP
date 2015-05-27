---------------------------------------
---------------------------------------
CMAP (The Connectivity Map)
---------------------------------------
---------------------------------------

---------------------------------------
The Connectivity Map: Using Gene-Expression Signatures to Connect Small Molecules, Genes, and Disease
Justin Lamb et al. Science 313, 1929 (2006); DOI: 10.1126/science.1132939
---------------------------------------

---------------------------------------
Source files downloaded from ftp://ftp.broad.mit.edu/pub/cmap/rankMatrix.txt.zip
---------------------------------------

---------------------------------------
---------------------------------------
Command line arguments for calculation
---------------------------------------
---------------------------------------


----------------------------
A.	source("path/cmap.R")
----------------------------

This command loads the R-file. The path has to be set if you load the connectivity
map from an extern folder. Leave the path away if the cmap-file is in the same folder
as you're working in. (Example: source("/media/bam23651_G/FinalVersion/cmap.R"))

---------------------------------------
B.	cmap(upRegFile, downRegFile, path)
---------------------------------------

The cmap is executed with three parameters.

upRegFile | downRegFile
	(required) Here the files are set, where the up and down regulated genes are
	saved. The files have to be ".grp"-files.

	upRegFile - The file with the up regulated genes. (Example: "readUp.grp")
	downRegFile - The file with the down regulated genes. (Example: "readDown.grp")

path
	(optional) If you work on Windows (normaly) the path has to be set. Or you're working
	on Linux or Mac then you can set the path to load your files from an extern folder.

	path - The path where the files can be found. (Example: "/media/bam23651_G/FinalVersion/")

----------------------------------------------
C. cmap.default(upRegFile, downRegFile, path)
----------------------------------------------

This function has the same use as in B. This is the default-operator for object-oriented 
programming in R.

---------------
D. cmap.init()
---------------

Here the rankmatrix is loaded. It is used for almost every calculation in this programm.
The function will be called automatically for cmap(...), cmap.default(...) and cmap.install(...).

------------------------------------------------
E. cmap.detResults(upRegFile, downRegFile, path)
------------------------------------------------

This command executes the calculation for the detailed results. It is loaded automatically
for cmap(...) and cmap.default(...).

upRegFile | downRegFile
	(required) Here the files are set, where the up and down regulated genes are
	saved. The files have to be ".grp"-files.

	upRegFile - The file with the up regulated genes. (Example: "readUp.grp")
	downRegFile - The file with the down regulated genes. (Example: "readDown.grp")

path
	(optional) If you work on Windows (normaly) the path has to be set. Or you're working
	on Linux or Mac then you can set the path to load your files from an extern folder.

	path - The path where the files can be found. (Example: "/media/bam23651_G/FinalVersion/")

------------------------------------------------------------
F. cmap.permResults(detRes, byNameORCellLine, path, install)
------------------------------------------------------------

This command executes the calculation for the permuted results. It is loaded automatically
for cmap(...) and cmap.default(...). It can only be calculated if you've calculated the 
detailed results before.

detRes | byNameORCellLine
	(required) Here the files are set, where the up and down regulated genes are
	saved. The files have to be ".grp"-files.

	detRes - The permuted results need the detailed results from E. for calculation.
	byNameORCellLine - Here you have to specify if you'd like to have the permuted results
	calculated by name or by name and cell line. There are only two choices.
	(Example: a.) "byName", b.) "CellLine")

path
	(optional) If you work on Windows (normaly) the path has to be set. Or you're working
	on Linux or Mac then you can set the path to load your files from an extern folder.

	path - The path where the files can be found. (Example: "/media/bam23651_G/FinalVersion/")

install
	(unusable) This parameter is set while installing the connectivity map. It has no use
	for a normal call. If you call the function with this parameter you won't observe a
	solution for the specificity.

	- install - Set if installation is on or off. There are only two choices.
	(Example: a.) "no", b.) "install")

--------------------------------------
G. cmap.install(path, LinuxORWindows)
--------------------------------------

With these command you'll install the cmap-files in a specific folder.

path
	(optional) If you work on Windows (normaly) the path has to be set. Or you're working
	on Linux or Mac then you can set the path to load your files from an extern folder and
	save the generated files in the new folder.

	path - The path where the files are installed. (Example: "/media/bam23651_G/FinalVersion/")

LinuxORWindows
	(optional) Here you can specify in which operating system you are working at. Linux is
	set as a default value.

	- LinuxORWindows - Sets the operating system. (Example: a.) "Linux", b.) "Windows")


----------------------------------------
----------------------------------------
Command line arguments for the solution
----------------------------------------
----------------------------------------


In this topic we will only see how to get the results after executing cmap(...). (When you
execute cmap.detResults you only have to save the output in an variable.
Example: matrix <- cmap.detResults)

detRes | nameRes | cellRes
	(required) With these commands you'll get your solutions.

	- detRes  - With detRes you'll get the detailed results. Possible calling methods:
				-	detRes[] => shows the whole matrix.
				-	detRes[i,] => shows the i-th row of the detailed results.
				-	detRes[i:j,] => shows the row from i to j of the detailed results.
				-	detRes["metformin",] => shows the first row where the instance metformin is listed.
				-	detRes[,i] => shows the i-th column of the detailed results.
				-	detRes[,i:j] => shows the column from i to j of the detailed results.
					(There are three columns: 1.) instance ID, 2.) score 3.) up regulation, 4.) down regulation)
				-	detRes[,"Score"] => shows the whole column with all scores of the detailed results.
					(other possibilities: "InstID", "Up", "Down")
				-	detRes[i, j] => get the j-th column of the i-th row.
				-	EXAMPLE:
					-	detRes[1:10, 2] => a vector with the scores of the first ten rows.

	If you'd like to have all results for a specific instance (e.g. "metformin", "geldanamycin", etc.) then copy and use this:
	x <- which(rownames(detRes) == "metformin")
	detRes[x,]

	- nameRes - With nameRes you'll get the permuted results calculated by name. Possible calling methods:
				-	nameRes[] => shows the whole matrix.
				-	nameRes[i,] => shows the i-th row of the permuted results.
				-	nameRes[i:j,] => shows the row from i to j of the permuted results.
				-	nameRes["metformin",] => shows the row where the instance metformin is listed.
				-	nameRes[,i] => shows the i-th column of the permuted results.
				-	nameRes[,i:j] => shows the column from i to j of the permuted results.
					(There are six columns: 1.) mean, 2.) n 3.) enrichment, 4.) p, 5.) specificity, 6.) % non-null)
				-	nameRes[,"enrichment"] => shows the whole column with all scores of the permuted results.
					(other possibilities: "mean", "n", "p", "specificity", "% non-null")
				-	nameRes[i, j] => get the j-th column of the i-th row.
				-	EXAMPLE:
					-	nameRes[1:10, 4] => a vector with the p values of the first ten rows.

	- cellRes - With nameRes you'll get the permuted results calculated by name and cell line. Possible calling methods:
				-	cellRes[] => shows the whole matrix.
				-	cellRes[i,] => shows the i-th row of the permuted results.
				-	cellRes[i:j,] => shows the row from i to j of the permuted results.
				-	cellRes["metformin - MCF7",] => shows the row where the instance metformin with its cell line is listed.
					(Here it is important to write " - " between the instance and the cell line!)
				-	cellRes[,i] => shows the i-th column of the permuted results.
				-	cellRes[,i:j] => shows the column from i to j of the permuted results.
					(There are six columns: 1.) mean, 2.) n 3.) enrichment, 4.) p, 5.) specificity, 6.) % non-null)
				-	cellRes[,"enrichment"] => shows the whole column with all scores of the permuted results.
					(other possibilities: "mean", "n", "p", "specificity", "% non-null")
				-	cellRes[i, j] => get the j-th column of the i-th row.
				-	EXAMPLE:
					-	cellRes[1:10, 4] => a vector with the p values of the first ten rows.


------------------------------
------------------------------
Files of the connectivity map
------------------------------
------------------------------


--------------
A. Folder
--------------

Four folders with data of the msig data base

msigdbUpByName | msigdbDownByName | msigdbUpByNameAndCell | msigdbDownByNameAndCell
	In these folders the data of the msig data base are saved.

	msigdbUpByName 			- Here the msig .grp-files for up regulation by name are saved.
	msigdbDownByName 		- Here the msig .grp-files for down regulation by name are saved.
	msigdbUpByNameAndCell	- Here the msig .grp-files for up regulation by name and cell line are saved.
	msigdbDownByNameAndCell	- Here the msig .grp-files for down regulation by name and cell line are saved.

--------------
B. .data-files
--------------

cellLineList.data | cmapData.data | instanceAndCellLineList.data | instanceCellNames.data | instanceList.data | instanceNames.data | instID.data
	These files hold the whole list of instances and cell lines. Without these files a calculation is not possible.

	cellLineList.data 				- The 6100 (current CMAP-Version) cell lines sorted by their instance ID.
	cmapData.data 					- Holds the whole cmap instance information. This file is needed for installation.
	instanceAndCellLineList.data 	- The 6100 (current CMAP-Version) instances with their cell lines sorted by their instance ID.
	instanceCellNames.data 			- The instances with their cell lines which are unique. Every entry is listed just once.
	instanceList.data 				- The 6100 (current CMAP-Version) instances sorted by their instance ID.
	instanceNames.data 				- The instances which are unique. Every entry is listed just once.
	instID.data 					- The 6100 (current CMAP-Version) instance ID's sorted by itself (instance ID).

--------------
C. .grp-files
--------------

readUp.grp | readDown.grp | alzheimerUp.grp | alzheimerDown.grp
	The first two files are just files for practicing and programming. They helped me to implement these programm.
	The second two files are from a scientific context. They were used for a scientific research.

	- readUp.grp 		- Holds two genes for up regulation.
	- readDown.grp 		- Holds two genes for down regulation.
	- alzheimerUp.grp 	- Holds genes for up regulation for a research for alzheimer.
	- alzheimerDown.grp - Holds genes for down regulation for a research for alzheimer.

--------------
D. .rda-files
--------------

rankmatrix.rda | pValue.rda | pValueByCell.rda | specificity.rda | specificityCell.rda
	These files save big arrays or matrices where the calculation would take too much time.

	- rankmatrix.rda 		- Holds the rankmatrix.
	- pValue.rda 			- Holds randList for the permuted results by name.
	- pValueByCell.rda 		- Holds randList for the permuted results by name and cell line.
	- specificity.rda 		- Holds specData for the permuted results by name.
	- specificityCell.rda 	- Holds specData for the permuted results by name and cell line.

--------------
E. .gmt-files
--------------

msigdb_up.gmt | msigdb_down.gmt
	These files can be downloaded from the cmap website. They hold the whole data of the msig data base query.

	- msigdb_up.gmt 	- Holds all gene sets from the msig data base for up regulation.
	- msigdb_down.gmt 	- Holds all gene sets from the msig data base for down regulation.

--------------
F. .txt-files
--------------

	rankMatrix.txt.sed1d | rankMatrix.txt.headn1

	- rankMatrix.txt.sed1d  - File of Anton Moll to install the rankmatrix.
	- rankMatrix.txt.headn1 - File of Anton Moll to install the rankmatrix.

--------------
G. .R-files
--------------

cmap.R
	The code of the connectivity map.


------------------------
------------------------
Examples from the shell
------------------------
------------------------


__________________________________________________________________________________

> detRes[1:10,]
                          InstID     Score        Up       Down
trimethylcolchicinic acid   2146 1.0000000 0.9411659 -0.8443657
diflorasone                 2142 0.9950989 0.9196697 -0.8571108
norfloxacin                 1406 0.9848694 0.9605978 -0.7979177
naproxen                    2533 0.9810742 0.8393843 -0.9123547
fluphenazine                2697 0.9806469 0.8720998 -0.8788763
benzethonium chloride       2508 0.9745394 0.8663106 -0.8737603
droperidol                  1290 0.9737854 0.8476866 -0.8910380
harmol                      1750 0.9704175 0.8633487 -0.8693623
LY-294002                   2696 0.9579260 0.8333707 -0.8770363
pentoxifylline              1444 0.9557142 0.8955706 -0.8108872

Detailed results for a query showing the first ten rows.
__________________________________________________________________________________

> nameRes["geldanamycin", ]
       mean           n  enrichment           p specificity  % non-null 
 -0.6481550  15.0000000  -0.6673224   0.0000000   0.2265193  93.3333333 

Permuted results for a query showing the result for the instance "geldanamycin".
__________________________________________________________________________________

> x <- which(cellRes[,"p"] <= 0.0005)
> cellRes[x,]
                          mean  n enrichment       p specificity % non-null
tanespimycin - PC3  -0.7527488 12 -0.8695082 0.00000 0.000000000  100.00000
monorden - MCF7     -0.7089552 12 -0.7537158 0.00000 0.000000000   91.66667
tanespimycin - MCF7 -0.7135461 36 -0.7366849 0.00000 0.008196721   94.44444
alvespimycin - MCF7 -0.7406855  7 -0.8519672 0.00001 0.000000000  100.00000
geldanamycin - MCF7 -0.6637222 10 -0.7262295 0.00003 0.022556391  100.00000
LY-294002 - MCF7    -0.3727243 34 -0.3605593 0.00015 0.368750000   58.82353

With this query you can get every instance with its cell line which its p value
is smaller or equaled 0.0005.
__________________________________________________________________________________


-------------------------
-------------------------
Design Decisions & Issues
-------------------------
-------------------------

Q. How to install?

A. For the installation you need to do only a few things:
	1.	Download the file "cmap_instances_02.xls" from the download page of the CMAP website ("instance inventory")
	2.	Open the file and mark the whole table (CTRL + A). And then copy the table (CTRL + C).
	3.	Create an empty file in your editor and paste the copied table (CTRL + V). Save the file as "cmapData.data" in your directory (CTRL + S).
	4.	Create four folders in your directory: "msigdbUpByName", "msigdbDownByName", "msigdbUpByNameAndCell", "msigdbDownByNameAndCell".
	5.	Don't forget to save the rankmatrix files: "rankMatrix.txt.headn1" and "rankMatrix.txt.sed1d".
	6.	Call the function for installation cmap.install(path, LinuxORWindows) where the path is your directory and LinuxORWindows is the operating system
		you are working at.

Q. In which order are all the objects and methods called?

A. Class Hierachy: the class hierarchy is found in the following file "hierarchy.pdf"