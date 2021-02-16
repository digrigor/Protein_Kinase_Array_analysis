library(data.table)
library(dplyr)


###USEFUL Functions####
uniprot_mapping <- function(ids) {
  uri <- 'http://www.uniprot.org/uniprot/?query='
  idStr <- paste(ids, collapse="+or+")
  format <- '&format=tab'
  fullUri <- paste0(uri,idStr,format)
  dat <- read.delim(fullUri, stringsAsFactors = F)
  dat
}


###Input masterfiles manipulation:


#User Defined input
names_mfile = "Protein Name_Masterfile.csv"
quant_mfile = "quantification_masterfile.csv"
uniprot_to_add_file = "ExtraProteinList.csv"

#User Defined output
network_files = paste0(dirname(quant_mfile),"../Network_files/")

#Input files
##uniprot names to ids
uniprot_mapping = unique(read.table("../Input/PhosphoDBs/uniprot-mapping.txt",stringsAsFactors = F, header = 1))

##phosphositeplus
pspl = data.frame(fread("../Input/PhosphoDBs/Kinase_Substrate_Dataset.csv"),stringsAsFactors = F)

##omni
omni = data.frame(fread("../Input/PhosphoDBs/Omnipath_enzyme_substrate_interactions_dataset.csv"),stringsAsFactors = F)

#Input files processing
metrics = read.csv(quant_mfile, stringsAsFactors = F)
uniprot_to_add = read.table(uniprot_to_add_file, stringsAsFactors = F)
row.names(uniprot_mapping)=make.names(uniprot_mapping$From,unique = T)
names = read.csv(names_mfile,stringsAsFactors = F)
names = names[,c(1:5)]
project = unlist(strsplit(basename(quant_mfile),"\\."))[1]

names$Membrane_Coordinate = gsub(" ","",names$Membrane_Coordinate)
metrics$Membrane.Coordinates = gsub(" ","",metrics$Membrane.Coordinates)

unip_to_add=uniprot_to_add$V1
#unip_to_add=c()

#Merge the protein names masterfile with the quantification masterfile tables.
nmetrics = merge(names,metrics,by.x="Membrane_Coordinate",by.y="Membrane.Coordinates",all.x = T)

#Create a new column (named "full_mod") that will be uniprot_ID-Modification (e.g. P23445-T354)
nmetrics$full_mod = ifelse(nmetrics$Modification!="",paste0(nmetrics$Uniprot_ID,"-",nmetrics$Modification),paste0(nmetrics$Uniprot_ID,"-"))
nmetrics$full_mod = gsub(" ","",nmetrics$full_mod)

#Subset the merged table and only get the columns we are interested in
#which are: Membrane_Coordinate, ArrayID.x, Uniprot_Name, Uniprot_ID, Modification, log2FC, full_mod
tka = nmetrics[,c(1:5,(ncol(nmetrics)-1),ncol(nmetrics))]

#Get the names of the proteins that do not have a phosphorylation site
no_kinases = c(tka[endsWith(tka$full_mod, "-"),]$Uniprot_ID, unip_to_add)


###Phosphorylation Databases manipulation:

##Omnipath kinase-substrate DB
#From the table get which molecule is the kinase and which is the substrate
omni$kinase = unlist(lapply(omni$name, function(x) unlist(strsplit(x, " "))[1]))
omni$substrate = unlist(lapply(omni$name, function(x) unlist(strsplit(x, " "))[4]))
omni$residue = paste0(omni$residue_type,omni$residue_offset)
omni$full_mod = paste0(omni$substrate,"-",omni$residue_type,omni$residue_offset)
omni$full_mod = gsub(" ","",omni$full_mod)

#Filter omni to only have phosphorylation/dephosphorylation interactions:
omni = subset(omni, modification %in% c("phosphorylation","dephosphorylation"))

#Subset the omni table to only get the columns important for us:
#which are:"SUID","kinase","modification","substrate","full_mod","source"
omni = omni[,c("SUID","kinase","modification","substrate","residue", "full_mod","sources")]




##PhosphositePlus kinase-substrate DB


pspl$kinase = pspl$KIN_ACC_ID
pspl$substrate = pspl$SUB_ACC_ID
pspl$residue = pspl$SUB_MOD_RSD
pspl$full_mod = paste0(pspl$substrate,"-",pspl$residue)
pspl$full_mod = gsub(" ","",pspl$full_mod)
 
#Filter PhosphositePlus to only have human molecules:
pspl = subset(pspl, KIN_ORGANISM == "human" & SUB_ORGANISM == "human")

#Subset the phosphosite table to only get the columns important for us:
#which are:"SUID","kinase","modification","substrate","full_mod","source"
pspl$modification = "phosphorylation"
pspl$SUID = sprintf("PSPL%05d",seq.int(nrow(pspl)))
pspl$sources = "PhosphositePlus"
pspl = pspl[,c("SUID","kinase","modification","substrate","residue", "full_mod","sources")]
colnames(pspl) <- colnames(omni)

##Merge the 2 databases together:
phosphodbs = rbind(omni,pspl)



###Create the network files
mynetwork = subset(phosphodbs,
                  (full_mod %in% tka$full_mod & kinase %in% unique(c(tka$Uniprot_ID,unip_to_add))) |
                  (substrate %in% unique(c(tka$Uniprot_ID,unip_to_add)) & kinase %in% unique(c(tka$Uniprot_ID,unip_to_add))))
mynetwork$kinase_genename = uniprot_mapping[mynetwork$kinase,]$To
mynetwork$substrate_genename = uniprot_mapping[mynetwork$substrate,]$To

##Make sure to separate multiple phosphorylation sites
nlist=c()
qlist=c()
for(i in c(1:nrow(mynetwork))){
  row = mynetwork[i,]
  kin = row$kinase
  logfcs = nmetrics[nmetrics$Uniprot_ID==kin,]$log2FC
  if(length(unique(logfcs))>1){
    print(row$kinase_genename)
    grlogs = nmetrics[nmetrics$Uniprot_ID==kin,c("log2FC","Modification")] %>%
      arrange(Modification) %>%
      group_by(log2FC) %>% 
      summarise(mods = paste(Modification, collapse = "/")) %>%
      data.frame(stringsAsFactors = F)
      for(i in c(1:nrow(grlogs))){
        grlog = grlogs[i,]
        row_alt = row
        row_alt$kinase = paste(c(kin,grlog$mods),collapse = "_")
        row_alt$kinase_genename = paste(c(row$kinase_genename,grlog$mods),collapse = "_")
        nlist = c(nlist,row_alt)
        qlist = c(qlist, row_alt$kinase, row_alt$kinase_genename, grlog$log2FC)
      }
  } else {
    nlist = c(nlist,row)
    qlist = c(qlist, row$kinase, row$kinase_genename, logfcs[1])
  }
}

mynetwork_kinases_fixed = as.data.frame(matrix(unlist(nlist),byrow = T, ncol = 9), stringsAsFactors = F)
colnames(mynetwork_kinases_fixed) = colnames(mynetwork) 

##Make sure to separate multiple phosphorylation sites
nlist=c()
for(i in c(1:nrow(mynetwork_kinases_fixed))){
  row = mynetwork_kinases_fixed[i,]
  kin = row$substrate
  logfcs = nmetrics[nmetrics$Uniprot_ID==kin,]$log2FC
  if(length(unique(logfcs))>1){
    print(row$substrate_genename)
    grlogs = nmetrics[nmetrics$Uniprot_ID==kin,c("log2FC","Modification")] %>%
      arrange(Modification) %>%
      group_by(log2FC) %>% 
      summarise(mods = paste(Modification, collapse = "/")) %>%
      data.frame(stringsAsFactors = F)
    if(row$residue %in% grlogs$mods){
      grlog = grlogs[grlogs$mods == row$residue,][1,]
      row_alt = row
      row_alt$substrate = paste(c(kin,grlog$mods),collapse = "_")
      row_alt$substrate_genename = paste(c(row$substrate_genename,grlog$mods),collapse = "_")
    } else {
      row_alt = row
    }
    nlist = c(nlist,row_alt)
    qlist = c(qlist, row_alt$substrate, row_alt$substrate_genename, grlog$log2FC)
  } else {
    row_alt = row
    nlist = c(nlist,row_alt)
    qlist = c(qlist, row$substrate, row$substrate_genename, logfcs[1])
  }
}

#Write final output
mynetwork_fixed = as.data.frame(matrix(unlist(nlist),byrow = T, ncol = 9), stringsAsFactors = F)
colnames(mynetwork_fixed) = colnames(mynetwork) 
quant_fixed = unique(as.data.frame(matrix(qlist,byrow = T, ncol=3),stringsAsFactors = F))

mynetwork_to_write = mynetwork_fixed[,c("kinase_genename","kinase","modification","substrate_genename","substrate","sources")]

write.csv(mynetwork_to_write, paste0(network_files,project,"_extra_network.txt"), quote = F, row.names = F)
write.csv(quant_fixed, paste0(network_files,project,"_extra_quantification.txt"), quote = F, row.names = F)
write.csv(phosphodbs, paste0("../Input/PhosphoDBs/phosphoDB.csv"), quote = F, row.names = F)
