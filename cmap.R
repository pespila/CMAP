#Attention: drugs = instances (cmap name)


#Defines a class as a function with 3 parameters.
#upRegfile: File with the upregulated genes; downRegFile: File with the downregulated genes; path: where your files are saved.
cmap <- function(upRegfile, downRegFile, path = "") UseMethod("cmap")

#Loads the rankmatrix if it's not loaded yet. Sets "rankmatrix" as a public-parameter.
#path: the directory
cmap.init <- function(path = "") {
        if (exists("rankmatrix") == FALSE) {
                #Get the file by path
                rankmatrixFile <- paste(path, "rankmatrix.rda", sep="")
                load(rankmatrixFile)
                rankmatrix <<- rankmatrix
        }
}

#This is the default setting of CMAP. You'll get the whole results (detailled, permuted) and some infos how you can call them.
#upRegfile: File with the upregulated genes; downRegFile: File with the downregulated genes; path: the directory
cmap.default <- function(upRegfile, downRegFile, path = "") {
        #With a time modul which gives you the seconds the calculation took.
        start <- proc.time()
        #get the results as matrices
        detRes <<- cmap.detResults(upRegfile, downRegFile, path) #detailed results
        nameRes <<- cmap.permResults(detRes, "byName", path) #permuted results by name
        cellRes <<- cmap.permResults(detRes, "CellLine", path) #permuted results by name and cell line

        stop <- proc.time()
        #calculate the running time
        time <- round(as.numeric(stop[3]-start[3]), digits=2)
        #Help and other widgets
        onScreen <- c("Your computing time was: ", time, " seconds.")
        cat(paste(strwrap(onScreen, 40)))
        cat("\n")

        getAccess <- c("01 Access: ", "To get access to your results type 'detRes', 'nameRes' or 'cellRes'.", "\n")
        cat(paste(strwrap(getAccess, 40)))
        cat("\n")

        helpWithAccess <- c("02 Help: ", "Type e.g. detRes[1:10,], nameRes[1:10,] or cellRes[1:10,] to see the first ten, best scored instances (in cellRes with the Cell Line).", "\n",
                "Or detRes['metformin',], nameRes['geldanamycin',] or cellRes['tanespimycin - MCF7'] to see the the solution for this instance.", "\n",
                "For a specific result in the detailed results use: 'InstID', 'Score', 'Up', 'Down'; in the permuted results: 'mean', 'n', 'enrichment', 'p', 'specificity', '% non-null'.", "\n",
                "e.g detRes[1:10, 'InstID'] (vecotr) or nameRes['geldanamycin', 'enrichment'] (value)", "\n")
        cat(paste(strwrap(helpWithAccess, 40)))
        cat("\n")
}

#This method calculates the detailled results of CMAP.
#upRegfile: File with the upregulated genes; downRegFile: File with the downregulated genes; path: where your files are saved.
cmap.detResults <- function(upRegfile, downRegFile, path = "") {
        cmap.init(path)
        #Add the path if it is set. If not -> the file name stays the same
        upRegPath <- paste(path, upRegfile, sep="")
        downRegPath <- paste(path, downRegFile, sep="")

        #This function is based on the Kolmogorov-Smirnoff-Statistic. It calculates the maximum distance between given distributions and the uniform distribution.
        detailedMaxValue <- function(v, n) {
                t <- length(v)
                j <- 1:t
                a <- max(j/t - v/n)
                b <- max(v/n - (j-1)/t)
                if (a > b) {k <- a}
                else {k <- -b}
                return (k) #returns the k_is => maximum distance between two distributions
        }

        #The detailed results function. It calculates all essential values and saves it in a matrix which is the return value.
        #rankmatrix: the loaded rankmatrix; upRegPath: File with the upregulated genes; downRegPath: File with the downregulated genes; path: where your files are saved.
        detResult <- function(rankmatrix, upRegPath, downRegPath, path = "") {
                #Get the number of instances and genes (n) from the rankmatrix
                instances <- ncol(rankmatrix)
                n <- nrow(rankmatrix)

                #Reads every line of a given file with its path. Theese files include the genes you're interested in.
                upReg <- as.character(readLines(upRegPath))
                downReg <- as.character(readLines(downRegPath))

                #Get the file by path
                IDFile <- paste(path, "instID.data", sep="")
                instanceFile <- paste(path, "instanceList.data", sep="")

                #Reads the instance ID, name of the instance, name of the cell line. The last two are connected to the rownames of the output-matrix.
                instID <- as.numeric(readLines(IDFile))
                instanceNames <- as.character(readLines(instanceFile))

                #Matrix with the scores of the genes your interested in.
                upRegMat <- rankmatrix[upReg, ]
                downRegMat <- rankmatrix[downReg, ]

                #Initialize lists
                scores <- c()
                ksUp <- c()
                ksDown <- c()

                #Calculation of the score.
                for (i in 1:instances) {
                        #Vectors for the calculation of the maximum distances in the detailedMaxValue function. Here the up regulated gene values are used.
                        v_up <- sort(upRegMat[,i])
                        k_up <- detailedMaxValue(v_up, n)

                        #Vectors for the calculation of the maximum distances in the detailedMaxValue function. Here the down regulated gene values are used.
                        v_down <- sort(downRegMat[,i])
                        k_down <- detailedMaxValue(v_down, n)

                        #Set s_i value by comparing the k values.
                        if (k_up >= 0 && k_down >= 0 || k_up < 0 && k_down < 0) {s_i <- 0}
                        else {s_i <- k_up-k_down}

                        #Save calculated values in the lists.
                        ksUp[i] <- k_up
                        ksDown[i] <- k_down
                        scores[i] <- s_i
                }

                #Get the maximum/minimum of all scores and implement the formular
                p <- max(scores)
                q <- min(scores)
                for (i in 1:instances) {
                        if (scores[i] >= 0) {scores[i] <- scores[i]/p}
                        else {scores[i] <- -(scores[i]/q)}
                }

                #Save lists in a matrix as columns. Give matrix the rownames and every column a colname. At last order the matrix by scores and ksUp.
                detResults <- matrix(c(instID, scores, ksUp, ksDown), instances)
                rownames(detResults) <- instanceNames
                colnames(detResults) <- c("InstID", "Score", "Up", "Down")
                detResults <- detResults[order(detResults[,2], detResults[,3], decreasing=TRUE),]

                return (detResults) #returns the matrix with the detailled results
        }
        return (detResult(rankmatrix, upRegPath, downRegPath, path)) #Return of the method. Returns the detailled results
}

#This is the function which calculates the permuted results of CMAP (both: by Name and by Name and Cell line).
#detRes: takes the detailled results (scores) for calculation; path: sets the path where files are saved. install: only for installation (don't needed for calculation)
cmap.permResults <- function (detRes, byNameORCellLine, path = "", install = "no") {
        #This function is based on the Kolmogorov-Smirnoff-Statistic. It calculates the maximum distance between given distributions and the uniform distribution.
        permutedMaxValue <- function(v, n) {
                t <- length(v)
                j <- 1:t
                a <- max(j/t - v/n)
                b <- max(v/n - (j-1)/t)
                if (a > b) {k <- a}
                else {k <- -b}
                return (k)
        }

        #This function returns the row numbers of a matrix where the rowname equals a given query
        getRowNumberByRowName <- function(mat, rowName) {
                rowNumbers <- which(rownames(mat) == rowName)
                return (rowNumbers) #returns a vector/list of indice
        }

        #The nonNull-function of CMAP returns a value. It proves how the scores are distributed in the detailed results. It returns the percentage how often the score is NOT NULL.
        #Example: 5 times above NULL, 4 times under NULL and 2 times Null. => 11 scores. 5 above and 4 under Null => take 5/11*100 = nonNull value.
        nonNull <- function(vec) {
                #proof where the given vector elements are greater or less Null
                lessNull <- which(vec < 0) #vector with numbers
                greaterNull <- which(vec > 0) #vector with numbers

                #Proof with if-clauses where the conditions are full filled and set your value
                if (length(lessNull) == 0 && length(greaterNull) == 0) {nonNull <- 0}
                if (length(lessNull) > length(greaterNull)) {nonNull <- length(lessNull)}
                else {nonNull <- length(greaterNull)}

                #get in percent
                nonNull <- nonNull*100/length(vec)

                return (nonNull) #returns a value
        }

        #For the permuted results by name and cell line you need other rownames. This function returns the detailed results as calculated but with
        #different rownames (instance name + cell line)
        addCellLine <- function(matrix, path = "") {
                #order the given matrix by instance ID
                orderByInstID <- matrix[order(matrix[,1]),]

                #Get files by path - Here the instance names and cell lines are read and the put together.
                IDFile <- paste(path, "cellLineList.data", sep="")
                instanceFile <- paste(path, "instanceList.data", sep="")
                cellLine <- as.character(readLines(IDFile))
                instanceNames <- as.character(readLines(instanceFile))
                rowNames <- paste(instanceNames, cellLine, sep=" - ")

                #Give the orderd matrix the new rownames
                rownames(orderByInstID) <- rowNames

                #Bring the matrix in the right order again.
                orderByScore <- orderByInstID[order(orderByInstID[,2], orderByInstID[,3], decreasing=TRUE),]

                return (orderByScore) #returns a matrix
        }

        permResult <- function (rankmatrix, detRes, byNameORCellLine, path, install) {
                #Get the number of instances from the rankmatrix
                instances <- ncol(rankmatrix)

                if(byNameORCellLine == "byName") {
                        #Get files by path
                        namesFile <- paste(path, "instanceNames.data", sep="")

                        #Reads all instances used in CMAP out of a file. Afterwards assign how much instances are used.
                        instanceNames <- as.character(readLines(namesFile))
                        instLen <- length(instanceNames)

                        #Get files by path
                        pFile <- paste(path, "pValue.rda", sep="")
                        specFile <- paste(path, "specificity.rda", sep="")

                        } else if (byNameORCellLine == "CellLine") {
                                #Get files by path
                                namesFile <- paste(path, "instanceCellNames.data", sep="")

                                #Reads all instances used in CMAP out of a file. Afterwards assign how much instances are used.
                                instanceNames <- as.character(readLines(namesFile))
                                instLen <- length(instanceNames)

                                #Get files by path
                                pValueFile <- "pValueByCell.rda"
                                specificityFile <- "specificityCell.rda"
                                pFile <- paste(path, pValueFile, sep="")
                                specFile <- paste(path, specificityFile, sep="")

                }

                #Initialize some lists.
                n <- c()
                p <- c()
                nonNull <- c()
                meanOfScores <- c()
                specificity <- c()
                enrichment <- c()

                #Load the p-test and specificity list
                load(pFile)
		if (install == "no") {load(specFile)}
		else if (install == "install") {specData <- matrix(0,instLen,100)} #only for the installation needed!!!

                #Main calculation
                for (i in 1:instLen) {
                        #Get the vector v for the permutedMaxValue function. Therefor you need the indices (rank) of the genes of the detailed results in descending order.
                        v <- getRowNumberByRowName(detRes, instanceNames[i])
                        instanceScores <- detRes[v, 2, drop = FALSE]

                        #The frequency of the genes which is set in the permuted results as n
                        n[i] <- length(v)

                        #The enrichment score is also calculated by the Kolmogorov-Smirnoff-Statistic.
                        enrichment[i] <- permutedMaxValue(v, instances) #uses permutedMaxValue function
                        
                        #Calculating mean and nonNull value             
                        meanOfScores[i] <- mean(instanceScores) #uses R function
                        nonNull[i] <- nonNull(instanceScores)       #uses nonNull function

                        #Calculates the p-value. Here randList of the pValue.rda/pValueByCell.rda loading is used.
                        #ATTENTION!!!!!!!! In the website based CMAP the p[i] for the if-state are set "---". Here 100 is used for 
                        #sorting the column correctly!!!
                        if (n[i] == 1 || meanOfScores[i] == 0 || nonNull[i] <= 50) {p[i] <- 100}
                        else {
                                p[i] <- sum(randList[[as.character(n[i])]] >= abs(enrichment[i]))/length(randList[[1]])
                        }

                        #Calculates the specificity. Here specData of specificity.rda/specificityCell.rda loading is used.
                        #The algorithm is based on the CMAP-Topic-File on the internet.
                        msigInstanceScore <- specData[i,]
                        if (enrichment[i] == 0) {
                                #If the calculated enrichment is Null set specificity Null
                                specificity[i] <- 0

                        }
                        else {
                                if (enrichment[i] < 0) {
                                        #If enrichment is less than Null get all indices of negative scores of the msigDB query with the given instance
                                        #("i-th" instance in the for-loop)
                                        negativeScores <- which(msigInstanceScore < 0)

                                        #get the scores which are less than Null
                                        getScores <- msigInstanceScore[negativeScores]

                                        #get indices where your computed enrichment is equaled or exceeded by the msig scores
                                        comparedScores <- which(enrichment[i] >= getScores)

                                        #the length of the compared scores devided by all negative msig scores is your specificity
                                        specificity[i] <- length(comparedScores)/length(negativeScores)

                                } else if (enrichment[i] > 0) {
                                        #Here the same computing is used except: now the programm is looking for all positive enrichment scores
                                        positiveScores <- which(msigInstanceScore > 0)
                                        getScores <- msigInstanceScore[positiveScores]
                                        comparedScores <- which(enrichment[i] >= getScores)
                                        specificity[i] <- 1-length(comparedScores)/length(positiveScores)

                                }

                        }
                }

                #Save lists in a matrix as columns. Give matrix the rownames and every column a colname. At last order the matrix by p-Value and enrichment.
                permutedResults <- matrix(c(meanOfScores, n, enrichment, p, specificity, nonNull),instLen)
                rownames(permutedResults) <- instanceNames
                colnames(permutedResults) <- c("mean", "n", "enrichment", "p", "specificity", "% non-null")
                permutedResults <- permutedResults[order(permutedResults[,4], permutedResults[,3], decreasing=FALSE),]

                return (permutedResults) #returns a matrix with the permuted results
        }

        #The most important query: result by name OR by name and cell line. With a wrong query the programm stops!!!
        if (byNameORCellLine == "byName") { #if the permuted results by name are wanted the method goes in this clause

                return (permResult(rankmatrix, detRes, byNameORCellLine, path, install))
        
        } else if (byNameORCellLine == "CellLine") { #if the permuted results by name and cell line are wanted the method goes in this clause
                
                return (permResult(rankmatrix, addCellLine(detRes, path), byNameORCellLine, path, install))
        
        }
        else {
                        stop("Only the parameters 'byName' OR 'CellLine' are possible.") #Here an error occurs.
        }
}

#Creates the files the programm needs.
#path: where your files are saved.
#LinuxORWindows: As it turns out that Linux and Windows use different seperators you have to specify in which operating system your work.
cmap.install <- function(path = "", LinuxORWindows = "Linux") {
        #By Anton Moll - Institue of Functional Genomics, path added by myself
        # this function requires generating two files derived from the original rankMatrix.txt (see comments below)
        readrankmatrix <- function(path = "") {
                # reading the whole matrix at once doesn't work, because read.table complains. so read header and matrix separately.
                # generate the two files from rankMatrix.txt, which can be downloaded from ftp://ftp.broad.mit.edu/pub/cmap/rankMatrix.txt.zip
                # execute these two commands:
                # # head -n1 rankMatrix.txt > rankMatrix.txt.headn1
                # # sed 1d rankMatrix.txt > rankMatrix.txt.sed1d
                head <- paste(path, "rankMatrix.txt.headn1", sep="")
                sed <- paste(path, "rankMatrix.txt.sed1d", sep="")
                rankmatrix.colnames <- read.table(head, sep="\t")      # read first line
                rankmatrix <- read.table(sed, sep="\t", row.names=1)   # read rest
                colnames(rankmatrix) <- rankmatrix.colnames[seq(2,6101)]
                rankmatrix <- as.matrix(rankmatrix)
                rankmatrix <<- rankmatrix
                saveFile <- paste(path, "rankmatrix.rda", sep="")
                save(rankmatrix, file = saveFile)
        }

        #--------------------------------------------------------------------------------

        #This function creates all the needed data-files.
        createDataFiles <- function(path = "") {
        		cmap.init(path)
                #get file by path
                cmapData <- paste(path, "cmapData.data", sep="")

                #read lines of all data and split it.
                linesOfAllData <- readLines(cmapData)
                splittedLines <- strsplit(linesOfAllData, split="\t", fixed=TRUE)

                #delete the first row -> in the xls file the first row.
                splittedLines <- splittedLines[-1]

                #initalize some lists
                instID <- c()
                cellLine <- c()
                instanceNames <- c()
                singleInstanceNames <- c()
                instanceAndCell <- c()

                #Save the values of the splitted data list and save it: 1 => Instance ID's (1.Column), 3 => instanceNames (3. Column), 7 => cell Lines (7. Column)
                for (i in 1:length(splittedLines)) {
                        instID[i] <- splittedLines[[i]][1]
                        cellLine[i] <- splittedLines[[i]][7]
                        instanceNames[i] <- splittedLines[[i]][3]
                }

                instances <- ncol(rankmatrix)

                x <- length(instID) #proof how long the list is. You could also use the cellLine list or the instanceName list.
                z <- instances+1 #Set a variable z to instances+1 (here: 6101)
                v <- z:x #Make a range from z to the real length

                #Cut off everything which hast nothing to do with instance names, cell lines or instance ID's
                instID <- instID[-v]
                cellLine <- cellLine[-v]
                instanceNames <- instanceNames[-v]

                #add the path
                instFile <- paste(path, "instID.data", sep="")
                cellLineFile <- paste(path, "cellLineList.data", sep="")
                instanceListFile <- paste(path, "instanceList.data", sep="")

                #Write the lists to files
                write(instID, file = instFile, ncolumns=1)
                write(cellLine, file = cellLineFile, ncolumns=1)
                write(instanceNames, file = instanceListFile, ncolumns=1)

                helpinstanceNames <- instanceNames #needed before instanceNames are deleted

                #This loop deletes instances which are listed several times. Goal is that every instance is listed only once!!!
                n <- 1
                while(length(instanceNames) > 0) {
                        singleInstanceNames[n] <- instanceNames[1]
                        x <- which(instanceNames == singleInstanceNames[n])
                        instanceNames <- instanceNames[-x]
                        n <- n + 1
                }

                #add the path
                instanceNamesFile <- paste(path, "instanceNames.data", sep="")

                #Save the instances names
                write(singleInstanceNames, file = instanceNamesFile, ncolumns=1)

                #a file is created which saves every gene with its cell line. It's needed for the file which saves only the unique genes with its cell line.
                a <- as.character(helpinstanceNames)
                b <- as.character(cellLine)
                instanceAndCellLineList <- paste(a, b, sep = " - ")

                #add path
                dacllFile <- paste(path, "instanceAndCellLineList.data", sep="")

                write(instanceAndCellLineList, file = dacllFile, ncolumns=1)

                #This loop creates a list where instance names with their cell line are saved
                n <- 1
                while(length(instanceAndCellLineList) > 0) {
                        instanceAndCell[n] <- instanceAndCellLineList[1]
                        x <- which(instanceAndCellLineList == instanceAndCell[n])
                        instanceAndCellLineList <- instanceAndCellLineList[-x]
                        n <- n + 1
                }

                #add path
                instanceCellNamesFile <- paste(path, "instanceCellNames.data", sep="")

                #Save the instance names with their cell line
                write(instanceAndCell, file = instanceCellNamesFile, ncolumns = 1)
        }

        #--------------------------------------------------------------------------------

        #The same function as you find in the permuted results ("getRowNumberByRowName")
        #The parameters changed but the function is still the same.
	numberOfInstancesByName <- function(vec, instanceName) {
                numberOfInstances <- which(vec == instanceName)
                return (numberOfInstances)
        }

        #returns how often a specific instance is listed in a vector
	getFrequency <- function (vec, instanceNames) {
                #initalize list
                instanceFrequ <- c()
                #proof the length for the loop
                len <- length(instanceNames)
                for (i in 1:len) {

                        #The length of this list gives you the frequency how often a instance is listed
                        lenOfInstances <- length(numberOfInstancesByName(vec, instanceNames[i]))

                        #save this number in your list
                        instanceFrequ[i] <- lenOfInstances
                }

                #With unique you get a list with different frequencies (same frequencies are deleted)
                instanceFrequ <- unique(instanceFrequ)

                return (instanceFrequ) #Returns a vector
        }

        #This function is based on the Kolmogorov-Smirnoff-Statistic. It calculates the maximum distance between given distributions and the uniform distribution.
        pMaxValue <- function(v, n) {
                t <- length(v)
                j <- 1:t
                a <- max(j/t - v/n)
                b <- max(v/n - (j-1)/t)
                if (a > b) {k <- a}
                else {k <- -b}
                return (k) #returns the k_is => Maximum Distance between two distributions
        }

        #The p-Value function. It returns a vector within enrichment scores of the instances by random ranks.
        pValue <- function (lens, maxLen) {
                #initalize a list
                pValue <- c()

                #The for-loop which has to count to d = 100.000 like the formula of CMAP says.
                #This sorting algorithm took less time than any other. For 100.000 calculations this is the best solution
                for (i in 1:100000) {
                        #Here random values are set
                        value <- sample.int(maxLen, lens)

                        #These random values are sorted by their indices
                        random <- sort.list(value, na.last = NA, method = "radix")

                        #And here the value list is sorted by the sorted indices => with this commands a random list is sorted.
                        value <- value[random]

                        #The enrichment is calculated with the pMaxValue function (=> the same as detailedMaxValue or permutedMaxValue)
                        pValue[i] <- pMaxValue(value, maxLen)
                }
                #For the algorithm only the absolute scores are needed
                pValue <- abs(pValue)

                return(pValue) #returns a vector
        }

        #This function saves the p-Values in an rda-File.
        savePValueFile <- function(byNameORCellLine, path = "") {
                #Initalize an empty list
                randList <- c()

                #get the number of instances which are saved in the rankmatrix
                instances <- ncol(rankmatrix)

                #proves if the solution should be calculated by name or by name and cell line and sets the right file for the next query
                if (byNameORCellLine == "byName") {
                        instanceNamesFile <- paste(path, "instanceNames.data", sep="")
						allInstancesFile <- paste(path, "instanceList.data", sep="")
                } else if (byNameORCellLine == "CellLine") {
                        instanceNamesFile <- paste(path, "instanceCellNames.data", sep="")
						allInstancesFile <- paste(path, "instanceAndCellLineList.data", sep="")
                }

                #read the instance names (with/or without their cell lines)
                instanceNames <- as.character(readLines(instanceNamesFile))
				allInstances <- as.character(readLines(allInstancesFile))

                #get the lengt of the instance names (with/or without their cell line)
                instLen <- length(instanceNames)

                #1. getFrequency: calculates how often the instance is saved in the rankmatrix
                #2. sort: sorts the list
                #3. [-1]: Without the first element (it is 1, so unusefull)
                var <- sort(getFrequency(allInstances, instanceNames))[-1]

                #in this loop every list of pValues (look at pValue-function) ist saved as a list in the randList (a List of random experiments)
                for (i in 1:length(var)) {
                        randList[[i]] <- pValue(var[i], instances)
                }

                #Every element of the randList (every list) gets a specific name => the number of Frequency
                names(randList) <- var

                #proves whether it's saved by Name or by Name and Cell Line.
                if (byNameORCellLine == "byName") {
                        pValueFile <- paste(path, "pValue.rda", sep="")
                } else if (byNameORCellLine == "CellLine") {
                        pValueFile <- paste(path, "pValueByCell.rda", sep="")
                }

                #saves the rda-file by Name/by Name and cell line with the right path
                save(randList, file = pValueFile)
        }

        #--------------------------------------------------------------------------------

        #calculates the specificity and saves it in an rda-file
        writeAndSaveSpec <- function(byNameORCellLine, path="", LinuxORWindows) {
                #sets the mode on installation
		install <- "install"

                #set the right seperator for your operating system
                if (LinuxORWindows == "Linux") {
                        sep <- "/"
                } else if (LinuxORWindows == "Windows") {
                        sep <- "\\"
                }

                #gets the gene sets with the upregulated genes
                msigUpFile <- paste(path, "msigdb_up.gmt", sep="")
                up <- as.character(readLines(msigUpFile))

                #gets the gene sets with the downregulated genes
                msigDownFile <- paste(path, "msigdb_down.gmt", sep="")
                down <- as.character(readLines(msigDownFile))

                #splits the list of the gene sets in a new list with lists
                u <- strsplit(up, split="\t", fixed=TRUE)
                d <- strsplit(down, split="\t", fixed=TRUE)
                lenghtOfMsigData <- length(u)

                #Initalize lists
                namesUp <- c()
                namesDown <- c()

                if (byNameORCellLine == "byName") {
                	upFolder <- "msigdbUpByName"
                	downFolder <- "msigdbDownByName"
					geneFile <- paste(path, "instanceNames.data", sep="")
                	numberOfGenes <- length(as.character(readLines(geneFile)))
                } else if (byNameORCellLine == "CellLine") {
                	upFolder <- "msigdbUpByNameAndCell"
                	downFolder <- "msigdbDownByNameAndCell"
					geneFile <- paste(path, "instanceCellNames.data", sep="")
                	numberOfGenes <- length(as.character(readLines(geneFile)))
                }

                #This loop just takes the list with the list of genes and saves the list of genes in an .grp-file. So later this grp-file is used for calculating.
                for (i in 1:lenghtOfMsigData) {
                        #Becuase the length(u)==length(d) it doesn't matter what you use. x <- u[[i]] just takes the first list in the list.
                        x <- u[[i]]
                        #The name of the probe of an upregulated set ist saved in a list
                        namesUp[i] <- x[1]
                        #deletes the first two elements of the list (name and empty space)
                        y <- x[-1:-2]
                        #the file name is created. ATTENTION: The right folder has to be created before running this loop!!!
                        upFile <- paste(path, upFolder, sep, "msigdbUp", i, ".grp", sep="")
                        #The files are saved in the folder
                        write(y, file = upFile, ncolumns=1)

                        #Here the same procedure is done for the down regulated gene sets
                        x <- d[[i]]
                        namesDown[i] <- x[1]
                        y <- x[-1:-2]                        
                        downFile <- paste(path, downFolder, sep, "msigdbDown", i, ".grp", sep="")
                        write(y, file = downFile, ncolumns=1)
                }

                #Initialize a matrix which saves all calculation of the MsigDB data (enrichment scores) for every gene (1309 genes with 312 MsigDB data)
                specData <- matrix(, numberOfGenes, lenghtOfMsigData)

                #in this loop the files are read (the saved files above) and the enrichment score (with detRes and permRes) is calculated as before and the
                #list with all scores of one query is saved in the matrix.
                for (i in 1:lenghtOfMsigData) {
                        #load the files without the path, because the path will be added by calling the detRes and permRes functions. Otherwise the path would
                        #be added two times and the programm wouldn't work.
                        loadUpFile <- paste(upFolder, sep, "msigdbUp", i, ".grp", sep="")
                        loadDownFile <- paste(downFolder, sep, "msigdbDown", i, ".grp", sep="")

                        detRes <- cmap.detResults(loadUpFile, loadDownFile, path)

                        if (byNameORCellLine == "byName") {
                          permRes <- cmap.permResults(detRes, "byName", path, install)
                        } else if (byNameORCellLine == "CellLine") {
                          permRes <- cmap.permResults(detRes, "CellLine", path, install)
                        }
                        specData[,i] <- permRes[,"enrichment"]
                }

                #The matrix is saved in an rda-file in the right folder.
                if (byNameORCellLine == "byName") {
                  saveFolder <- paste(path,"specificity.rda",sep="")
                } else if (byNameORCellLine == "CellLine") {
                  saveFolder <- paste(path,"specificityCell.rda", sep="")
                }
                save(specData, file = saveFolder)
        }

        #--------------------------------------------------------------------------------

        #Calls the functions with the right parameters
		print("The installation has started.")
        readrankmatrix(path) #Starts the function
		print("Rankmatrix rda-file installed.")

        createDataFiles(path) #Starts the function
		print("Data Files installed")

        savePValueFile("byName", path) #Starts the function
		print("p-Value rda-file 'by Name' installed.")

        savePValueFile("CellLine", path) #Starts the function
		print("p-Value rda-file 'by Name and Cell Line' installed.")

        writeAndSaveSpec("byName", path, LinuxORWindows) #Starts the function
		print("specificity rda-file 'by Name' completed.")

        writeAndSaveSpec("CellLine", path, LinuxORWindows) #Starts the function
		print("specificity rda-file 'by Name and Cell Line' installed.")
}