#' @noRd
deprintize<-function(f){
  return(function(...) {capture.output(w<-f(...));return(w);});
}

#' Run CCC methods
#'
#' Run cell-cell communication methods
#'
#' @param name_mat Name of scRNA-seq data
#' @param Lcell Name of ligand(sender) cell
#' @param Rcell Name of receptor(receiver) cell
#' @param mat_ori Original expression matrix
#' @param label_ori Original cell label
#' @param list_methods List of cell-cell communication methods
#' @param path_result Path of result
#'
#' @return
#' @import reticulate
#' @export
#'
#' @examples
run_CCC_methods <- function(name_mat, Lcell, Rcell, mat_ori, label_ori, list_methods, path_result){
  list_methods <- tolower(list_methods)
  path_result_CCC <- paste(path_result,'CCCmethods',sep='//')
  if(!dir.exists(path_result_CCC)){dir.create(path_result_CCC)}
  ####################################################################################
  ## 11 method based on R
  ####################################################################################
  # M1 CellCall
  if("cellcall" %in% list_methods){
    print("CellCall,start")
    deprintize(run_CellCall)(mat_ori, label_ori, path_result_CCC, name_mat)
    print("CellCall,end")
  }

  # M2 CellChat
  if("cellchat" %in% list_methods){
    print("CellChat,start")
    deprintize(run_CellChat)(mat_ori, label_ori, path_result_CCC, name_mat, Lcell, Rcell)
    print("CellChat,end")
  }
  # M3 CytoTalk
  if("cytotalk" %in% list_methods){
    print("CytoTalk,start")
    deprintize(run_Cytotalk)(mat_ori, label_ori, path_result_CCC, name_mat, Lcell, Rcell)
    print("CytoTalk,end")
  }

  # M4 iCellNet
  if("icellnet" %in% list_methods){
    print("iCellNet,start")
    deprintize(run_iCellNet)(mat_ori, label_ori, path_result_CCC, name_mat, Lcell, Rcell)
    print("iCellNet,end")
  }

  # M5 iTALK
  if("italk" %in% list_methods){
    print("iTALK,start")
    deprintize(run_iTalk)(mat_ori, label_ori, path_result_CCC, name_mat)
    print("iTALK,end")
  }

  # M6 NicheNet
  if("nichenet" %in% list_methods){
    print("NicheNet,start")
    deprintize(run_NicheNet)(mat_ori, label_ori, path_result_CCC, name_mat, Lcell, Rcell)
    print("NicheNet,end")
  }

  # M7 scMLnet
  if("scmlnet" %in% list_methods){
    print("scMLnet,start")
    deprintize(run_scMLnet)(mat_ori, label_ori, path_result_CCC, name_mat, Lcell, Rcell)
    print("scMLnet,end")
  }

  # M8 SingleCellSignalR
  if("singlecellsignalr" %in% list_methods){
    print("SingleCellSignalR,start")
    deprintize(run_SingleCellSignalR)(mat_ori, label_ori, path_result_CCC, name_mat, Lcell, Rcell)
    print("SingleCellSignalR,end")
  }

  # M9 Zhou
  if("zhou" %in% list_methods){
    print("Zhou,start")
    deprintize(run_Zhou)(mat_ori, label_ori, path_result_CCC, name_mat, Lcell, Rcell)
    print("Zhou,end")
  }

  # M10 Skelly
  if("skelly" %in% list_methods){
    print("Skelly,start")
    deprintize(run_Skelly)(mat_ori, label_ori, path_result_CCC, name_mat, Lcell, Rcell)
    print("Skelly,end")
  }

  # M11 Kumar
  if("kumar" %in% list_methods){
    print("Kumar,start")
    deprintize(run_Kumar)(mat_ori, label_ori, path_result_CCC, name_mat, Lcell, Rcell)
    print("Kumar,end")
  }

  ####################################################################################
  ## 3 method based on python
  ####################################################################################
  # M12 NATMI
  if("natmi" %in% list_methods){
    print("NATMI,start")
    natmi <- import_from_path('NATMI_py', path=system.file('python', package='RobustCCC'), convert=TRUE)
    path_NATMI_info = system.file('depend_cell_cell_communication_methods','NATMI',package='RobustCCC')
    deprintize(natmi$run_NATMI)(path_NATMI_info, mat_ori, label_ori, path_result_CCC, name_mat)
    print("NATMI,end")
  }

  # M13 scConnect
  if("scconnect" %in% list_methods){
    print("scConnect,start")
    scconnect <- import_from_path('scconnect_py', path=system.file('python', package='RobustCCC'), convert=TRUE)
    deprintize(scconnect$run_scConnect)(mat_ori, label_ori, path_result_CCC, name_mat)
    print("scConnect,end")
  }

  # M14 CellPhoneDB
  if("cellphonedb" %in% list_methods){
    print("CellPhoneDB,start")
    cellphonedb <- import_from_path('CellPhoneDB_py', path=system.file('python', package='RobustCCC'), convert=TRUE)
    path_cpdb = system.file('depend_cell_cell_communication_methods','CellPhoneDB',package='RobustCCC')
    deprintize(cellphonedb$run_cellphoneDB)(mat_ori, label_ori, path_result_CCC, path_cpdb, name_mat)
    print("CellPhoneDB,end")

  }
}


#' Run CCC methods in simulated data
#'
#' Load simulated data and run cell-cell communication methods in simulated data
#'
#' @param list_simulated_type List of simulated type
#' @param name_mat_CTP Name of scRNA-seq data for cell type pairs
#' @param Lcell Name of ligand(sender) cell
#' @param Rcell Name of receptor(receiver) cell
#' @param list_methods List of cell-cell communication methods
#' @param path_data Path of data
#' @param path_result Path of result
#'
#' @return
#' @export
#'
#' @examples
run_CCC_methods_simulated <-
  function(list_simulated_type, name_mat_CTP, Lcell, Rcell, list_methods, path_data, path_result){
  num_of_rep = 3
  for(simulated_type in list_simulated_type){
    if(simulated_type=="simuReplicate"){list_rate = c(95,90,85)}else{list_rate = c(5,10,15)}
    for(rate in list_rate){
      for(rep_i in 0:(num_of_rep-1)){

        # load simulated mat
        if(simulated_type %in% c('simuReplicate', 'GaussianNoise', 'dropout', 'selfNoised', 'ligRecPermu', 'expPermuted')){
          name_mat_CTP_simulate = paste(name_mat_CTP,simulated_type,rate,rep_i,sep='_')
        }else{name_mat_CTP_simulate =name_mat_CTP}
        mat_simulate = read.table(paste(path_data,name_mat_CTP_simulate,sep = '//'),sep=',',header = T,row.names = 1)
        mat_simulate = mat_simulate[!duplicated(rownames(mat_simulate)),]

        # change name_mat_CTP_simulate
        name_mat_CTP_simulate = paste(name_mat_CTP,simulated_type,rate,rep_i,sep='_')

        # load simulated label
        if(simulated_type %in% c('simuReplicate', 'cellTypePermu')){
          name_label_CTP_simulate = paste(name_label_CTP,simulated_type,rate,rep_i,sep='_')
        }else{name_label_CTP_simulate=name_label_CTP}
        label_simulate = read.table(paste(path_data,name_label_CTP_simulate,sep='//'),sep=',',header = TRUE, row.names = 1)

        # run cell-cell communication
        run_CCC_methods(name_mat_CTP_simulate, Lcell, Rcell, mat_simulate, label_simulate, list_methods, path_result)
      }
    }
  }
}

#' @noRd
run_CellCall <- function(mat_ori, label_ori, path_result_CCC, name_mat){

  path_write_curr = paste(path_result_CCC,'CellCall',sep='//')
  if(!dir.exists(path_write_curr)){dir.create(path_write_curr)}
  path_write_curr_error = paste(path_write_curr,'ERROR',sep='//')
  if(!dir.exists(path_write_curr_error)){dir.create(path_write_curr_error)}

  name_save = name_mat
  if(name_save %in% list.files(path_write_curr)){return(0)}
  if(name_save %in% list.files(path_write_curr_error)){return(0)}

  try_idx = 1
  try_idx = try({
    label_ori['cellType'] = label_ori['cell_type']

    cellName_old = colnames(mat_ori)
    cellName_ner = sapply(cellName_old, function(x) paste(x,gsub('_','',label_ori[x,'cellType']),sep='-'))
    in.content = mat_ori
    colnames(in.content) = cellName_ner


    mt <- cellcall::CreateNichConObject(data=in.content, min.feature = 0,
                                        names.field = 2,
                                        names.delim = "-",
                                        source = "TPM",
                                        Org = "Homo sapiens",
                                        project = "Microenvironment")# Homo sapiens, Mus musculus
    mt <- suppressMessages(cellcall::TransCommuProfile(object = mt,use.type="mean",Org = "Homo sapiens"))

    score = mt@data[["expr_l_r_log2_scale"]]
    colName_old = colnames(score)
    score['LR'] = rownames(score)
    score = score[,c('LR',colName_old)]

    if (length(rownames(score))>0){
      write.table(score, paste(path_write_curr,name_save,sep='//'), row.names=FALSE,col.names=TRUE,sep='\t',quote = FALSE)
    }
    else{write.table(score, paste(path_write_curr_error,name_save,sep='//'), row.names=FALSE,col.names=TRUE,sep='\t',quote = FALSE)}
    suppressWarnings(rm(mat_ori,label_ori, in.content, mt, score));gc()
  })
  if('try-error' %in% class(try_idx)){cat(try_idx, file = paste(path_write_curr_error,name_save,sep='//'), sep = "\n")}
  suppressWarnings(rm(mat_ori,label_ori, in.content, mt, score));gc()
}

#' @noRd
run_CellChat <- function(mat_ori, label_ori, path_result_CCC, name_mat, Lcell, Rcell){
  path_write_curr = paste(path_result_CCC,'CellChat',sep='//')
  if(!dir.exists(path_write_curr)){dir.create(path_write_curr)}
  path_write_curr_error = paste(path_write_curr,'ERROR',sep='//')
  if(!dir.exists(path_write_curr_error)){dir.create(path_write_curr_error)}

  name_save = name_mat
  if(name_save %in% list.files(path_write_curr)){return(0)}
  if(name_save %in% list.files(path_write_curr_error)){return(0)}

  try_idx = 1
  try_idx = try({
    # create cellchat object
    cellchat <- CellChat::createCellChat(object = as.matrix(mat_ori), meta = label_ori, group.by = "cell_type")
    cellchat@DB <- CellChat::CellChatDB.human # CellChatDB.human, CellChatDB.mouse

    # preprocessing
    cellchat <- CellChat::subsetData(cellchat)
    cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
    cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)

    # communication network
    cellchat <- CellChat::computeCommunProb(cellchat)
    net.prob = cellchat@net$prob
    net.pval = cellchat@net$pval

    LRnameList = dimnames(net.pval)[[3]]
    table_all  = data.frame()
    for(i in 1:length(LRnameList)){
      prob = net.prob[,,i][Lcell,Rcell]
      pval = net.pval[,,i][Lcell,Rcell]
      if(pval<1){
        table_temp = data.frame(LRname=LRnameList[i],prob=prob,pval=pval)
        table_all = rbind(table_all,table_temp)
      }
    }
    table_all = table_all[order(table_all[,3],table_all[,2]),]
    write.table(table_all, paste(path_write_curr,name_save,sep='//'), row.names=FALSE,col.names=TRUE,sep=',',quote = FALSE)
    suppressWarnings(rm(mat_ori,label_ori,cellchat,net.prob,net.pval,LRnameList,table_all));gc()
  })
  if('try-error' %in% class(try_idx)){cat(try_idx, file = paste(path_write_curr_error,name_save,sep='//'), sep = "\n")}
  suppressWarnings(rm(mat_ori,label_ori,cellchat,net.prob,net.pval,LRnameList,table_all));gc()
}

#' @noRd
run_Cytotalk <- function(mat_ori, label_ori, path_result_CCC, name_mat, Lcell, Rcell){

  path_write_temp = paste(path_result_CCC,'CytoTalk',sep='//')
  if(!dir.exists(path_write_temp)){dir.create(path_write_temp)}
  path_write_curr = paste(path_write_temp,name_mat,sep='//')
  if(!dir.exists(path_write_curr)){dir.create(path_write_curr)}

  if('Crosstalk_TypATypB.txt' %in% list.files(path_write_curr)){return(0)}

  try_idx = 1
  try_idx = try({

    label_ori['cellType'] = label_ori['cell_type']
    label_ori_type_a = label_ori[label_ori['cell_type']==Lcell,]
    mat_type_a = mat_ori[,rownames(label_ori_type_a)]
    label_ori_type_b = label_ori[label_ori['cell_type']==Rcell,]
    mat_type_b = mat_ori[,rownames(label_ori_type_b)]

    type_a = Lcell
    type_b = Rcell
    dir_out= path_write_curr
    proteins=CytoTalk::pcg_human
    ligands=CytoTalk::ligands_human
    cutoff_a=0.1;cutoff_b=0.1

    dir_out <- suppressWarnings(normalizePath(dir_out))

    # make sure output directory exists
    if (!dir.exists(dir_out)) {dir.create(dir_out)}

    preprocess(proteins, mat_type_a, mat_type_b, cutoff_a, cutoff_b, dir_out)
    compute_non_self_talk(ligands, mat_type_a, mat_type_b, dir_out)
    status <- compute_pem(mat_type_a, mat_type_b,type_a,type_b, dir_out)
    compute_crosstalk(type_a, type_b, dir_out)
    suppressWarnings(rm(proteins,ligands));gc()
  })
  if('try-error' %in% class(try_idx)){cat(try_idx, file = paste(path_write_curr,'ERROR',sep='//'), sep = "\n");return(0)}
  suppressWarnings(rm(proteins,ligands));gc()
}

#' @noRd
run_iCellNet <- function(mat_ori, label_ori, path_result_CCC, name_mat, Lcell, Rcell){

  iCellNet_db <- system.file('depend_cell_cell_communication_methods/iCellNet','ICELLNETdb.tsv',package='RobustCCC')
  db=as.data.frame(read.csv(iCellNet_db, sep="\t",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""))
  db.name.couple=icellnet::name.lr.couple(db, type="Family")

  path_write_curr = paste(path_result_CCC,'ICELLNET',sep='//')
  if(!dir.exists(path_write_curr)){dir.create(path_write_curr)}
  path_write_curr_error = paste(path_write_curr,'ERROR',sep='//')
  if(!dir.exists(path_write_curr_error)){dir.create(path_write_curr_error)}

  name_save = name_mat
  if(name_save %in% list.files(path_write_curr)){return(0)}
  if(name_save %in% list.files(path_write_curr_error)){return(0)}

  try_idx = 1
  try_idx = try({
    # 2 - Retrieve gene expression matrix
    data = mat_ori
    target = label_ori
    target['Cluster'] = label_ori['cell_type']
    target$Cell=rownames(target)

    average.manual=matrix(ncol=length(unique(target$Cluster)), nrow=length(rownames(data)))
    colnames(average.manual)=unique(target$Cluster)
    rownames(average.manual)=rownames(data)
    for (cell in unique(target$Cluster)){
      cells.clust=target$Cell[which(target$Cluster==cell)]
      average.manual[,cell]=apply(data[,which(colnames(data)%in%cells.clust)], 1, mean)
    }
    average.clean=average.manual

    # 3 - Apply icellnet pipeline on cluster of interest
    ## cc->PC
    data.icell=as.data.frame(icellnet::gene.scaling(as.data.frame(average.clean), n=1, db=db))
    PC.data=as.data.frame(data.icell[,c(Rcell, "Symbol")], row.names = rownames(data.icell))
    PC.target=data.frame("Class"=c(Rcell), "ID"= c(Rcell), "Cell_type"=c(Rcell))
    rownames(PC.target)=c(Rcell)
    my.selection=c(Rcell)

    score.computation.1= icellnet::icellnet.score(direction="out", PC.data=PC.data,
                                                  CC.data= as.data.frame(data.icell[,Lcell], row.names = rownames(data.icell)),
                                                  PC.target = PC.target, PC=my.selection, CC.type = "RNAseq",
                                                  PC.type = "RNAseq",  db = db)
    lr1=score.computation.1[[2]]
    lr1 = sort(lr1[(lr1[,1]>0)&(!is.na(lr1[,1])),],decreasing = TRUE)
    lr1_table = data.frame(lr=names(lr1),score=lr1)

    write.table(lr1_table, paste(path_write_curr,name_save,sep='//'), row.names=FALSE,col.names=TRUE,sep=',',quote = FALSE)
    suppressWarnings(rm(mat_ori,data,label_ori,target,average.manual,average.clean,data.icell,PC.data,PC.target,score.computation.1,lr1,lr1_table));gc()
  })
  if('try-error' %in% class(try_idx)){cat(try_idx, file = paste(path_write_curr_error,name_save,sep='//'), sep = "\n")}
  suppressWarnings(rm(mat_ori,data,label_ori,target,average.manual,average.clean,data.icell,PC.data,PC.target,score.computation.1,lr1,lr1_table));gc()
}

#' @noRd
run_iTalk <- function(mat_ori, label_ori, path_result_CCC, name_mat){

  path_write_curr = paste(path_result_CCC,'iTALK',sep='//')
  if(!dir.exists(path_write_curr)){dir.create(path_write_curr)}
  path_write_curr_error = paste(path_write_curr,'ERROR',sep='//')
  if(!dir.exists(path_write_curr_error)){dir.create(path_write_curr_error)}

  name_save = name_mat
  if(name_save %in% list.files(path_write_curr)){return(0)}
  if(name_save %in% list.files(path_write_curr_error)){return(0)}

  try_idx = 1
  try_idx = try({

    mat_ori = data.frame(t(mat_ori))
    mat_ori['cell_type'] = label_ori['cell_type']
    highly_exprs_genes<-iTALK::rawParse(mat_ori,top_genes=50,stats='mean')
    comm_list<-c('growth factor','other','cytokine','checkpoint')
    res<-NULL
    for(comm_type in comm_list){
      res_cat<-iTALK::FindLR(highly_exprs_genes,datatype='mean count',comm_type=comm_type)
      res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
      res<-rbind(res,res_cat)
    }

    write.table(res, paste(path_write_curr,name_save,sep='//'), row.names=FALSE,col.names=TRUE,sep=',',quote = FALSE)
    suppressWarnings(rm(mat_ori,label_ori,res));gc()
  })
  if('try-error' %in% class(try_idx)){cat(try_idx, file = paste(path_write_curr_error,name_save,sep='//'), sep = "\n")}
  suppressWarnings(rm(mat_ori,label_ori,res));gc()
}

#' @noRd
run_Kumar <- function(mat_ori, label_ori, path_result_CCC, name_mat, Lcell, Rcell){

  path_LR <- system.file('depend_cell_cell_communication_methods/Kumar_Skelly_Zhou','PairsLigRec_simple_human.txt',package='RobustCCC')
  table_LR = read.table(path_LR,sep='\t',header = T)

  path_write_curr = paste(path_result_CCC,'Kumar',sep='//')
  if(!dir.exists(path_write_curr)){dir.create(path_write_curr)}
  path_write_curr_error = paste(path_write_curr,'ERROR',sep='//')
  if(!dir.exists(path_write_curr_error)){dir.create(path_write_curr_error)}

  name_save = name_mat
  if(name_save %in% list.files(path_write_curr)){return(0)}
  if(name_save %in% list.files(path_write_curr_error)){return(0)}

  try_idx = 1
  try_idx = try({
    label_ori['cellType'] = label_ori['cell_type']

    table_LR_temp = table_LR
    table_LR_temp = table_LR_temp[(table_LR_temp$Ligand %in% rownames(mat_ori))
                                  &(table_LR_temp$Receptor %in% rownames(mat_ori)),]

    mat_L = mat_ori[table_LR_temp[table_LR_temp$Ligand %in% rownames(mat_ori),'Ligand'],
                    rownames(label_ori[label_ori['cell_type']==Lcell,])]
    mat_R = mat_ori[table_LR_temp[table_LR_temp$Ligand %in% rownames(mat_ori),'Receptor'],
                    rownames(label_ori[label_ori['cell_type']==Rcell,])]
    num_Lcells_Expressed = rowSums(mat_L)
    num_Rcells_Expressed = rowSums(mat_R)
    table_LR_temp[,c('exp_L','exp_R','product_LR')] = list(num_Lcells_Expressed,num_Rcells_Expressed,num_Lcells_Expressed*num_Rcells_Expressed)



    write.table(table_LR_temp, paste(path_write_curr,name_save,sep='//'), row.names=FALSE,col.names=TRUE,sep=',',quote = FALSE)
    suppressWarnings(rm(mat_ori,label_ori, mat_L, mat_R, table_LR_temp, num_Lcells_Expressed, num_Rcells_Expressed, table_LR));gc()
  })
  if('try-error' %in% class(try_idx)){cat(try_idx, file = paste(path_write_curr_error,name_save,sep='//'), sep = "\n")}
  suppressWarnings(rm(mat_ori,label_ori, mat_L, mat_R, table_LR_temp, num_Lcells_Expressed, num_Rcells_Expressed, table_LR));gc()

}

#' @import Seurat
run_NicheNet <- function(mat_ori, label_ori, path_result_CCC, name_mat, Lcell, Rcell){

  if(!file.exists(paste(system.file('depend_cell_cell_communication_methods/NicheNet',package='RobustCCC'),'ligand_target_matrix.rds',sep='//'))){
    ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
    # ligand_target_matrix = readRDS(url("https://media.githubusercontent.com/media/chenxing-zhang/RobustCCC/main/inst/depend_cell_cell_communication_methods/NicheNet/ligand_target_matrix.rds"))
    saveRDS(ligand_target_matrix,paste(system.file('depend_cell_cell_communication_methods/NicheNet',package='RobustCCC'),'ligand_target_matrix.rds',sep='//'))
  }
  else if(file.size(paste(system.file('depend_cell_cell_communication_methods/NicheNet',package='RobustCCC'),'ligand_target_matrix.rds',sep='//'))<100000000){
    ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
    # ligand_target_matrix = readRDS(url("https://media.githubusercontent.com/media/chenxing-zhang/RobustCCC/main/inst/depend_cell_cell_communication_methods/NicheNet/ligand_target_matrix.rds"))
    saveRDS(ligand_target_matrix,paste(system.file('depend_cell_cell_communication_methods/NicheNet',package='RobustCCC'),'ligand_target_matrix.rds',sep='//'))
  }
  else{
    ligand_target_matrix <- readRDS(system.file('depend_cell_cell_communication_methods/NicheNet','ligand_target_matrix.rds',package='RobustCCC'))
  }

  if(!file.exists(paste(system.file('depend_cell_cell_communication_methods/NicheNet',package='RobustCCC'),'lr_network.rds',sep='//'))){
    lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
    saveRDS(lr_network,paste(system.file('depend_cell_cell_communication_methods/NicheNet',package='RobustCCC'),'lr_network.rds',sep='//'))
  }
  else{
    lr_network <- readRDS(system.file('depend_cell_cell_communication_methods/NicheNet','lr_network.rds',package='RobustCCC'))
  }


  if(!file.exists(paste(system.file('depend_cell_cell_communication_methods/NicheNet',package='RobustCCC'),'weighted_networks.rds',sep='//'))){
    weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
    saveRDS(weighted_networks,paste(system.file('depend_cell_cell_communication_methods/NicheNet',package='RobustCCC'),'weighted_networks.rds',sep='//'))
  }
  else{
    weighted_networks <- readRDS(system.file('depend_cell_cell_communication_methods/NicheNet','weighted_networks.rds',package='RobustCCC'))
  }


  path_write_curr = paste(path_result_CCC,'NicheNet',sep='//')
  if(!dir.exists(path_write_curr)){dir.create(path_write_curr)}
  path_write_curr_error = paste(path_write_curr,'ERROR',sep='//')
  if(!dir.exists(path_write_curr_error)){dir.create(path_write_curr_error)}

  name_save = name_mat
  if(name_save %in% list.files(path_write_curr)){return(0)}
  if(name_save %in% list.files(path_write_curr_error)){return(0)}

  try_idx = 1
  try_idx = try({
    mat_t = t(mat_ori)
    label_ori['celltype'] = label_ori['cell_type']

    # create markers
    scRNA <- CreateSeuratObject(counts = mat_ori, meta.data=label_ori)
    scRNA@meta.data['celltype'] = scRNA@meta.data['cell_type']
    Idents(scRNA) = scRNA@meta.data['cell_type']
    markers <- FindMarkers(scRNA, ident.1 = Rcell, ident.2 = Lcell)
    # if('try-error' %in% class(markers)){return(0)}

    # Step 1: Define expressed genes in sender and receiver cell populations
    sender_ids = rownames(label_ori[label_ori$'cell_type'==Lcell,])
    receiver_ids = rownames(label_ori[label_ori$'cell_type'==Rcell,])
    if(length(sender_ids)==1){expressed_genes_sender = mat_t[sender_ids,]%>% .[. >= 1] %>% names()}
    else{expressed_genes_sender = mat_t[sender_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 1] %>% names()}
    if(length(receiver_ids)==1){expressed_genes_receiver = mat_t[receiver_ids,]%>% .[. >= 1] %>% names()}
    else{expressed_genes_receiver = mat_t[receiver_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 1] %>% names()}

    # Step 2: Define the gene set of interest and a background of genes
    geneset_oi = rownames(markers[((markers$p_val<0.05)&(markers$avg_log2FC>0)),])
    background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

    # Step 3: Define a set of potential ligands
    ligands = lr_network %>% pull(from) %>% unique()
    expressed_ligands = intersect(ligands,expressed_genes_sender)
    receptors = lr_network %>% pull(to) %>% unique()
    expressed_receptors = intersect(receptors,expressed_genes_receiver)
    lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)
    potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

    # Step 4: Perform NicheNet ligand activity analysis on the gene set of interest
    ligand_activities = nichenetr::predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
    best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

    # Follow-up analysis 1: Ligand-receptor network inference for top-ranked ligands
    # get the ligand-receptor network of the top-ranked ligands
    lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
    best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
    # get the weights of the ligand-receptor interactions as used in the NicheNet model
    lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
    lr_network_top_df =lr_network_top_df[order(lr_network_top_df$weight,decreasing=TRUE),]

    write.table(lr_network_top_df, paste(path_write_curr,name_save,sep='//'), row.names=FALSE,col.names=TRUE,sep=',',quote = FALSE)

    suppressWarnings(rm(mat_ori,label_ori,mat_t,scRNA,lr_network_top_df,markers));gc()
  })
  if('try-error' %in% class(try_idx)){cat(try_idx, file = paste(path_write_curr_error,name_save,sep='//'), sep = "\n")}
  suppressWarnings(rm(mat_ori,label_ori,mat_t,scRNA,lr_network_top_df,markers));gc()

}

#' @import Seurat
run_scMLnet <- function(mat_ori, label_ori, path_result_CCC, name_mat, Lcell, Rcell){
  pval=1;logfc=0;cores=1
  LigRecLib <- system.file('depend_cell_cell_communication_methods/scMLnet','LigRec.txt',package='RobustCCC')
  TFTarLib <- system.file('depend_cell_cell_communication_methods/scMLnet','TFTargetGene.txt',package='RobustCCC')
  RecTFLib <- system.file('depend_cell_cell_communication_methods/scMLnet','RecTF.txt',package='RobustCCC')


  path_write_curr = paste(path_result_CCC,'scMLnet',sep='//')
  if(!dir.exists(path_write_curr)){dir.create(path_write_curr)}
  path_write_curr_error = paste(path_write_curr,'ERROR',sep='//')
  if(!dir.exists(path_write_curr_error)){dir.create(path_write_curr_error)}

  name_save = name_mat
  if(name_save %in% list.files(path_write_curr)){return(0)}
  if(name_save %in% list.files(path_write_curr_error)){return(0)}

  try_idx = 1
  try_idx = try({
    mat_ori<- as(as.matrix(mat_ori),"dgCMatrix")
    label_ori['Barcode'] = rownames(label_ori)
    label_ori['Cluster'] = label_ori['cell_type']

    netList <- RunMLnet(mat_ori, label_ori, Lcell, Rcell,
                        pval, logfc,
                        LigRecLib, TFTarLib, RecTFLib,cores)

    LigRecTable = data.frame(t(data.frame(strsplit(netList$LigRec,"_"))))
    colnames(LigRecTable) = c('Lig','Rec')
    rownames(LigRecTable) = 1:dim(LigRecTable)[1]
    LigGenePavlue = netList$LigGenePavlue[LigRecTable[,'Lig'],]
    RecGenePavlue = netList$RecGenePavlue[LigRecTable[,'Rec'],]
    LigRecTable[,c('LogFC_Lig','Pva_Lig')] = LigGenePavlue[,c('LogFC','Pval')]
    LigRecTable[,c('LogFC_Rec','Pva_Rec')] = RecGenePavlue[,c('LogFC','Pval')]

    write.table(LigRecTable, paste(path_write_curr,name_save,sep='//'), row.names=FALSE,col.names=TRUE,sep=',',quote = FALSE)
    suppressWarnings(rm(mat_ori,label_ori,netList,LigGenePavlue, RecGenePavlue, LigRecTable));gc()
  })
  if('try-error' %in% class(try_idx)){cat(try_idx, file = paste(path_write_curr_error,name_save,sep='//'), sep = "\n")}
  suppressWarnings(rm(mat_ori,label_ori,netList,LigGenePavlue, RecGenePavlue, LigRecTable));gc()
}

#' @noRd
run_SingleCellSignalR <- function(mat_ori, label_ori, path_result_CCC, name_mat, Lcell, Rcell){

  library(SingleCellSignalR)
  path_write_curr = paste(path_result_CCC,'SingleCellSignalR',sep='//')
  if(!dir.exists(path_write_curr)){dir.create(path_write_curr)}
  path_write_curr_error = paste(path_write_curr,'ERROR',sep='//')
  if(!dir.exists(path_write_curr_error)){dir.create(path_write_curr_error)}

  name_save = name_mat
  if(name_save %in% list.files(path_write_curr)){return(0)}
  if(name_save %in% list.files(path_write_curr_error)){return(0)}

  try_idx = 1
  try_idx = try({
    label_ori['celltype'] = label_ori['cell_type']
    label_ori[label_ori['cell_type']==Lcell,'celltype']=1
    label_ori[label_ori['cell_type']==Rcell,'celltype']=2

    signal <- SingleCellSignalR::cell_signaling(data = mat_ori, genes = rownames(mat_ori), cluster = label_ori$celltype, write = FALSE,s.score=0,tol = 1)

    signal.12 = signal$`cluster 1-cluster 2`

    write.table(signal.12, paste(path_write_curr,name_save,sep='//'), row.names=FALSE,col.names=TRUE,sep=',',quote = FALSE)
    suppressWarnings(rm(mat_ori,label_ori,clust.ana,signal,signal.12));gc()
  })
  if('try-error' %in% class(try_idx)){cat(try_idx, file = paste(path_write_curr_error,name_save,sep='//'), sep = "\n")}
  suppressWarnings(rm(mat_ori,label_ori,clust.ana,signal,signal.12));gc()
}

#' @noRd
run_Skelly <- function(mat_ori, label_ori, path_result_CCC, name_mat, Lcell, Rcell){

  path_LR <- system.file('depend_cell_cell_communication_methods/Kumar_Skelly_Zhou','PairsLigRec_simple_human.txt',package='RobustCCC')
  table_LR = read.table(path_LR,sep='\t',header = T)

  path_write_curr = paste(path_result_CCC,'Skelly',sep='//')
  if(!dir.exists(path_write_curr)){dir.create(path_write_curr)}
  path_write_curr_error = paste(path_write_curr,'ERROR',sep='//')
  if(!dir.exists(path_write_curr_error)){dir.create(path_write_curr_error)}

  name_save = name_mat
  if(name_save %in% list.files(path_write_curr)){return(0)}
  if(name_save %in% list.files(path_write_curr_error)){return(0)}

  try_idx = 1
  try_idx = try({
    label_ori['cellType'] = label_ori['cell_type']
    table_LR_temp = table_LR
    table_LR_temp = table_LR_temp[(table_LR_temp$Ligand %in% rownames(mat_ori))
                                  &(table_LR_temp$Receptor %in% rownames(mat_ori)),]

    mat_L = mat_ori[table_LR_temp[table_LR_temp$Ligand %in% rownames(mat_ori),'Ligand'],
                    rownames(label_ori[label_ori['cell_type']==Lcell,])]
    mat_R = mat_ori[table_LR_temp[table_LR_temp$Ligand %in% rownames(mat_ori),'Receptor'],
                    rownames(label_ori[label_ori['cell_type']==Rcell,])]
    num_Lcells_Expressed = rowSums(mat_L>0)
    num_Rcells_Expressed = rowSums(mat_R>0)
    table_LR_temp[,c('num_Lcells_Expressed','rate_Lcells_Expressed')] = list(num_Lcells_Expressed,num_Lcells_Expressed/dim(mat_L)[2])
    table_LR_temp[,c('num_Rcells_Expressed','rate_Rcells_Expressed')] = list(num_Rcells_Expressed,num_Rcells_Expressed/dim(mat_R)[2])

    write.table(table_LR_temp, paste(path_write_curr,name_save,sep='//'), row.names=FALSE,col.names=TRUE,sep=',',quote = FALSE)
    suppressWarnings(rm(mat_ori,label_ori, mat_L, mat_R, table_LR_temp, num_Lcells_Expressed, num_Rcells_Expressed, table_LR));gc()
  })
  if('try-error' %in% class(try_idx)){cat(try_idx, file = paste(path_write_curr_error,name_save,sep='//'), sep = "\n")}
  suppressWarnings(rm(mat_ori,label_ori, mat_L, mat_R, table_LR_temp, num_Lcells_Expressed, num_Rcells_Expressed, table_LR));gc()

}

#' @import limma
run_Zhou<- function(mat_ori, label_ori, path_result_CCC, name_mat, Lcell, Rcell){

  path_LR <- system.file('depend_cell_cell_communication_methods/Kumar_Skelly_Zhou','PairsLigRec_simple_human.txt',package='RobustCCC')
  table_LR = read.table(path_LR,sep='\t',header = T)

  path_write_curr = paste(path_result_CCC,'Zhou',sep='//')
  if(!dir.exists(path_write_curr)){dir.create(path_write_curr)}
  path_write_curr_error = paste(path_write_curr,'ERROR',sep='//')
  if(!dir.exists(path_write_curr_error)){dir.create(path_write_curr_error)}

  name_save = name_mat
  if(name_save %in% list.files(path_write_curr)){return(0)}
  if(name_save %in% list.files(path_write_curr_error)){return(0)}

  try_idx = 1
  try_idx = try({
    label_ori['cellType'] = label_ori['cell_type']

    label_ori['Lcell'] = label_ori['cell_type']
    label_ori[label_ori['cell_type']==Lcell,'Lcell']=1
    label_ori[label_ori['cell_type']==Rcell,'Lcell']=0
    label_ori['Rcell'] = label_ori['cell_type']
    label_ori[label_ori['cell_type']==Lcell,'Rcell']=0
    label_ori[label_ori['cell_type']==Rcell,'Rcell']=1

    design = cbind(Lcell=as.numeric(label_ori[,'Lcell']),Rcell=as.numeric(label_ori[,'Rcell']))
    rownames(design) = rownames(label_ori)
    fit = lmFit(mat_ori,design);
    cont.matrix = makeContrasts(contrasts = c('Lcell-Rcell','Rcell-Lcell'),levels=design)
    fit2 = contrasts.fit(fit,cont.matrix);
    fit2 = eBayes(fit2)
    DEtable_L2R=topTable(fit2,coef = 'Lcell-Rcell', n = dim(mat_ori)[1]);
    DEtable_R2L=topTable(fit2,coef = 'Rcell-Lcell', n = dim(mat_ori)[1]);

    table_LR_temp = table_LR
    table_LR_temp = table_LR_temp[(table_LR_temp$Ligand %in% rownames(DEtable_L2R))
                                  &(table_LR_temp$Receptor %in% rownames(DEtable_R2L)),]

    table_LR_temp[,c('LogFC_Lig','Pva_Lig')] = DEtable_L2R[table_LR_temp$Ligand,c('logFC','adj.P.Val')]
    table_LR_temp[,c('LogFC_Rec','Pva_Rec')] = DEtable_R2L[table_LR_temp$Receptor,c('logFC','adj.P.Val')]

    write.table(table_LR_temp, paste(path_write_curr,name_save,sep='//'), row.names=FALSE,col.names=TRUE,sep=',',quote = FALSE)
    suppressWarnings(rm(mat_ori,label_ori,design, fit, cont.matrix, fit2, DEtable_L2R, DEtable_R2L, table_LR));gc()
  })
  if('try-error' %in% class(try_idx)){cat(try_idx, file = paste(path_write_curr_error,name_save,sep='//'), sep = "\n")}
  suppressWarnings(rm(mat_ori,label_ori,design, fit, cont.matrix, fit2, DEtable_L2R, DEtable_R2L, table_LR));gc()
}

