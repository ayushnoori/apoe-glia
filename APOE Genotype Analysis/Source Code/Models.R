## *****************************************************************************
## Model 1
## exp ~ APOE + CERAD/Braak

run_model1 = function(sf, mDataG, SF, E, annot, ds, celltype, useVoom = F) {
  
  # Regression Analysis
  # factors for creating design matrix 
  # SF: CERAD or Braak stages
  # E: APOE2 (including 22, 23), APOE3 (including 33), APOE4 (including 44, 34, 24)

  # remove samples with NA meta data
  sel = (!is.na(SF) & !is.na(E))
  SF = SF[sel]
  E = E[sel]
  mDataG = mDataG[,sel]
  cat(sum(!sel), "out of", length(sel), "samples were removed because of lack of clinical records\n")
  cat("Final sample size: ", sum(sel), "\n")
  
  # design matrix
  design.sfE = model.matrix(~ E + SF + 0 )
  colnames(design.sfE) = gsub("EE","E", colnames(design.sfE))
  
  
  # limma model
  fit.sfE = lmFit(mDataG, design.sfE)
  fit.c.sfE =  contrasts.fit(fit.sfE, makeContrasts(C1="E4-E3", C2="E2-E3", levels=design.sfE))
  if(useVoom){ 
    fit.eb.sfE = eBayes(fit.c.sfE, robust=TRUE)
  }else{
    fit.eb.sfE = eBayes(fit.c.sfE, trend = T)
  }
  
  # annotation
  gen.annot = merge(data.table(EnsemblID = rownames(mDataG)), annot, by = "EnsemblID", all.x = T, sort = F)
  fit.eb.sfE$genes = gen.annot
  
  # call topTable
  tT1.sfE = topTable(fit.eb.sfE, sort.by="none", number=Inf, coef="C1", confint=TRUE)
  tT2.sfE = topTable(fit.eb.sfE, sort.by="none", number=Inf, coef="C2", confint=TRUE)
  
  colnames(tT1.sfE)[3:10] = paste0("E4vsE3_", colnames(tT1.sfE)[3:10])
  colnames(tT2.sfE)[3:10] = paste0("E2vsE3_", colnames(tT2.sfE)[3:10])

  setDT(tT1.sfE)
  setDT(tT2.sfE)

  tT.sfE = tT1.sfE[tT2.sfE, on=c("GeneSymbol", "EnsemblID")]
  
  return(tT.sfE)
}

## *****************************************************************************
## Model 1.5
## exp ~ E4Num + E2Num + CERAD: 

run_model1.5 = function(sf, mDataG, SF, E4Num, E2Num, annot, ds, celltype, useVoom = F) {
  
  # regression Analysis
  # factors for creating design matrix 
  # E4Num: the number of copy of e4 (numeric)
  # E2Num: the number of copy of e2 (numeric)
  
  # remove samples with NA meta data
  sel = (!is.na(SF) & !is.na(E4Num) & !is.na(E2Num) )
  SF = SF[sel]
  E4Num = E4Num[sel]
  E2Num = E2Num[sel]
  mDataG = mDataG[,sel]
  cat(sum(!sel), "out of", length(sel), "samples were removed because of lack of clinical records\n")
  cat("Final sample size: ", sum(sel), "\n")
  
  # design matrix
  design.sfE = model.matrix(~ E4Num + E2Num + SF )
  
  # limma model
  fit.sfE = lmFit(mDataG, design.sfE)
  if(useVoom){ 
    fit.eb.sfE = eBayes(fit.sfE, robust=TRUE)
  }else{
    fit.eb.sfE = eBayes(fit.sfE, trend = T)
  }
  
  # annotation
  gen.annot = merge(data.table(EnsemblID = rownames(mDataG)), annot, by = "EnsemblID", all.x = T, sort = F)
  fit.eb.sfE$genes = gen.annot
  
  return(fit.eb.sfE)
}


## *****************************************************************************
## Model 1.5.2
## exp ~ E4Num + E2Num: for C0 subjects only

run_model1.5.2 = function(mDataG, E4Num, E2Num, annot, ds, celltype, useVoom = F) {
  
  # Regression Analysis
  # factors for creating design matrix 
  # E4Num: the number of copy of e4 (numeric)
  # E2Num: the number of copy of e2 (numeric)
  
  # remove samples with NA meta data
  sel = (!is.na(E4Num) & !is.na(E2Num) )
  E4Num = E4Num[sel]
  E2Num = E2Num[sel]
  mDataG = mDataG[,sel]
  cat(sum(!sel), "out of", length(sel), "samples were removed because of lack of clinical records\n")
  cat("Final sample size: ", sum(sel), "\n")
  
  # design matrix
  design.sfE = model.matrix(~ E4Num + E2Num)
  
  # limma model
  fit.sfE = lmFit(mDataG, design.sfE)
  if(useVoom){ 
    fit.eb.sfE = eBayes(fit.sfE, robust=TRUE)
  }else{
    fit.eb.sfE = eBayes(fit.sfE, trend = T)
  }
  
  # annotation
  gen.annot = merge(data.table(EnsemblID = rownames(mDataG)), annot, by = "EnsemblID", all.x = T, sort = F)
  fit.eb.sfE$genes = gen.annot
  
  return(fit.eb.sfE)
}


## *****************************************************************************
## Model 2
## baseline: C0E4-C0E3, C0E2-C0E3 | B1E4-B1E3, B1E2-B1E3
run_model2 = function(sf, mDataG, SF, E, annot, ds, celltype, useVoom = F) {
  
  # remove samples with NA meta data
  sel = (!is.na(SF) & !is.na(E))
  SF = SF[sel]
  E = E[sel]
  mDataG = mDataG[,sel]
  cat(sum(!sel), "out of", length(sel), "samples were removed because of lack of clinical records\n")
  cat("Final sample size: ", sum(sel), "\n")
  
  # factors for creating design matrix: SF.E
  # SF: CERAD or Braak stages
  # E: APOE2 (including 22, 23), APOE3 (including 33), APOE4 (including 44, 34, 24)
  if(sf=="CERAD") {
    f3=factor(paste0("C",SF,E))
  } else {
    f3=factor(paste0("B",SF,E))
  }
  
  # design matrix
  design.EBase = model.matrix(~f3 + 0)
  colnames(design.EBase) = gsub("f3","", colnames(design.EBase))
  
  # limma model
  if(sf=="CERAD") {
    fit.EBase = lmFit(mDataG, design.EBase)
    fit.c.EBase =  contrasts.fit(fit.EBase, makeContrasts(C1="C0E4-C0E3", C2="C0E2-C0E3", levels=design.EBase))
  } else {
    fit.EBase = lmFit(mDataG, design.EBase)
    fit.c.EBase =  contrasts.fit(fit.EBase, makeContrasts(C1="B1E4-B1E3", C2="B1E2-B1E3", levels=design.EBase))
  }
  
  if(useVoom){ 
    fit.eb.EBase = eBayes(fit.c.EBase, robust=TRUE)
  }else{
    fit.eb.EBase = eBayes(fit.c.EBase, trend = T)
  }
  
  # annotation
  gen.annot = merge(data.table(EnsemblID = rownames(mDataG)), annot, by = "EnsemblID", all.x = T, sort = F)
  fit.eb.EBase$genes = gen.annot
  
  # call topTable
  tT1.EBase = topTable(fit.eb.EBase, sort.by="none", number=Inf, coef="C1", confint=TRUE)
  tT2.EBase = topTable(fit.eb.EBase, sort.by="none", number=Inf, coef="C2", confint=TRUE)
  
  colnames(tT1.EBase)[3:10] = paste0("E4vsE3_", colnames(tT1.EBase)[3:10])
  colnames(tT2.EBase)[3:10] = paste0("E2vsE3_", colnames(tT2.EBase)[3:10])
  
  setDT(tT1.EBase)
  setDT(tT2.EBase)
  
  tT.EBase = tT1.EBase[tT2.EBase, on=c("GeneSymbol", "EnsemblID")]

  return(tT.EBase)
}


## *****************************************************************************
## Model 3
## C3E4-C3E3, C3E2-C3E3 | B3E4-B3E3, B3E2-B3E3
run_model3 = function(sf, mDataG, SF, E, annot, ds, celltype, useVoom = F) {
  
  # remove samples with NA meta data
  sel = (!is.na(SF) & !is.na(E))
  SF = SF[sel]
  E = E[sel]
  mDataG = mDataG[,sel]
  cat(sum(!sel), "out of", length(sel), "samples were removed because of lack of clinical records\n")
  cat("Final sample size: ", sum(sel), "\n")
  
  # factors for creating design matrix: SF.E
  # SF: CERAD or Braak stages
  # E: APOE2 (including 22, 23), APOE3 (including 33), APOE4 (including 44, 34, 24)
  if(sf=="CERAD") {
    f3=factor(paste0("C",SF,E))
  } else {
    f3=factor(paste0("B",SF,E))
  }
  
  design.EBase <- model.matrix(~f3 + 0)
  colnames(design.EBase) <- gsub("f3","", colnames(design.EBase))
  
  # limma model
  if(sf=="CERAD") {
    fit.EBase = lmFit(mDataG, design.EBase)
    fit.c.EBase =  contrasts.fit(fit.EBase, makeContrasts(C1="C3E4-C3E3", C2="C3E2-C3E3", levels=design.EBase))
  } else {
    fit.EBase = lmFit(mDataG, design.EBase)
    fit.c.EBase =  contrasts.fit(fit.EBase, makeContrasts(C1="B3E4-B3E3", C2="B3E2-B3E3", levels=design.EBase))
  }
  
  if(useVoom){ 
    fit.eb.EBase = eBayes(fit.c.EBase, robust=TRUE)
  }else{
    fit.eb.EBase = eBayes(fit.c.EBase, trend = T)
  }
  
  # annotation
  gen.annot = merge(data.table(EnsemblID = rownames(mDataG)), annot, by = "EnsemblID", all.x = T, sort = F)
  fit.eb.EBase$genes = gen.annot
  
  # call topTable
  tT1.EBase = topTable(fit.eb.EBase, sort.by="none", number=Inf, coef="C1", confint=TRUE)
  tT2.EBase = topTable(fit.eb.EBase, sort.by="none", number=Inf, coef="C2", confint=TRUE)
  colnames(tT1.EBase)[3:10] = paste0("E4vsE3_", colnames(tT1.EBase)[3:10])
  colnames(tT2.EBase)[3:10] = paste0("E2vsE3_", colnames(tT2.EBase)[3:10])
  
  setDT(tT1.EBase)
  setDT(tT2.EBase)
  
  tT.EBase = tT1.EBase[tT2.EBase, on=c("GeneSymbol", "EnsemblID")]
  
  return(tT.EBase)
}


## *****************************************************************************
## Model 4
## C23E4-C23E3, C23E2-C23E3
run_model4 = function(sf, mDataG, SF, E, annot, ds, celltype, useVoom = F) {
  
  # remove samples with NA meta data
  sel = (!is.na(SF) & !is.na(E))
  SF = SF[sel]
  E = E[sel]
  mDataG = mDataG[,sel]
  cat(sum(!sel), "out of", length(sel), "samples were removed because of lack of clinical records\n")
  cat("Final sample size: ", sum(sel), "\n")
  
  # factors for creating design matrix: SF.E
  # SF: CERAD or Braak stages
  # E: APOE2 (including 22, 23), APOE3 (including 33), APOE4 (including 44, 34, 24)
  SF = gsub("2|3", "23", SF)
  f3 = factor(paste0("C",SF,E))
  
  design.EBase = model.matrix(~f3 + 0)
  colnames(design.EBase) = gsub("f3","", colnames(design.EBase))
  
  # limma model
  if(sf=="CERAD") {
    fit.EBase = lmFit(mDataG, design.EBase)
    fit.c.EBase =  contrasts.fit(fit.EBase, makeContrasts(C1="C23E4-C23E3", C2="C23E2-C23E3", levels=design.EBase))
  } else {
    stop("only for CERAD")
  }
  
  if(useVoom){ 
    fit.eb.EBase = eBayes(fit.c.EBase, robust=TRUE)
  }else{
    fit.eb.EBase = eBayes(fit.c.EBase, trend = T)
  }
  
  # annotation
  gen.annot = merge(data.table(EnsemblID = rownames(mDataG)), annot, by = "EnsemblID", all.x = T, sort = F)
  fit.eb.EBase$genes = gen.annot
  
  # call topTable
  tT1.EBase = topTable(fit.eb.EBase, sort.by="none", number=Inf, coef="C1", confint=TRUE)
  tT2.EBase = topTable(fit.eb.EBase, sort.by="none", number=Inf, coef="C2", confint=TRUE)
  colnames(tT1.EBase)[3:10] = paste0("E4vsE3_", colnames(tT1.EBase)[3:10])
  colnames(tT2.EBase)[3:10] = paste0("E2vsE3_", colnames(tT2.EBase)[3:10])
  
  setDT(tT1.EBase)
  setDT(tT2.EBase)
  
  tT.EBase = tT1.EBase[tT2.EBase, on=c("GeneSymbol", "EnsemblID")]
  
  return(tT.EBase)
}


## *****************************************************************************
## Model 5
## CXE4-CXE3, CXE2-CXE3, where X could be 1 or 2
run_model5 = function(sf, mDataG, SF, E, annot, ds, celltype, useVoom = F) {
  
  # remove samples with NA meta data
  sel = (!is.na(SF) & !is.na(E))
  SF = SF[sel]
  E = E[sel]
  mDataG = mDataG[,sel]
  cat(sum(!sel), "out of", length(sel), "samples were removed because of lack of clinical records\n")
  cat("Final sample size: ", sum(sel), "\n")
  
  # factors for creating design matrix: SF.E
  # SF: CERAD or Braak stages
  # E: APOE2 (including 22, 23), APOE3 (including 33), APOE4 (including 44, 34, 24)
  f3 = factor(paste0("C",SF,E))
  
  design.EBase = model.matrix(~f3 + 0)
  colnames(design.EBase) = gsub("f3","", colnames(design.EBase))
  
  # limma model
  if(sf=="CERAD") {
    fit.EBase = lmFit(mDataG, design.EBase)
    fit.c.EBase =  contrasts.fit(fit.EBase, makeContrasts(C1="C1E4-C1E3", C2="C1E2-C1E3", 
                                                          C3="C2E4-C2E3", C4="C2E2-C2E3",levels=design.EBase))
  } else {
    stop("only for CERAD")
  }
  
  if(useVoom){ 
    fit.eb.EBase = eBayes(fit.c.EBase, robust=TRUE)
  }else{
    fit.eb.EBase = eBayes(fit.c.EBase, trend = T)
  }
  
  # annotation
  gen.annot = merge(data.table(EnsemblID = rownames(mDataG)), annot, by = "EnsemblID", all.x = T, sort = F)
  fit.eb.EBase$genes = gen.annot
  
  # call topTable
  tT1.EBase = topTable(fit.eb.EBase, sort.by="none", number=Inf, coef="C1", confint=TRUE)
  tT2.EBase = topTable(fit.eb.EBase, sort.by="none", number=Inf, coef="C2", confint=TRUE)
  tT3.EBase = topTable(fit.eb.EBase, sort.by="none", number=Inf, coef="C3", confint=TRUE)
  tT4.EBase = topTable(fit.eb.EBase, sort.by="none", number=Inf, coef="C4", confint=TRUE)
  colnames(tT1.EBase)[3:10] = paste0("C1E4vsC1E3_", colnames(tT1.EBase)[3:10])
  colnames(tT2.EBase)[3:10] = paste0("C1E2vsC1E3_", colnames(tT2.EBase)[3:10])
  colnames(tT3.EBase)[3:10] = paste0("C2E4vsC2E3_", colnames(tT3.EBase)[3:10])
  colnames(tT4.EBase)[3:10] = paste0("C2E2vsC2E3_", colnames(tT4.EBase)[3:10])
  
  setDT(tT1.EBase)
  setDT(tT2.EBase)
  setDT(tT3.EBase)
  setDT(tT4.EBase)
  
  tT.EBase = tT1.EBase[tT2.EBase, on=c("GeneSymbol", "EnsemblID")]
  tT.EBase = tT.EBase[tT3.EBase, on=c("GeneSymbol", "EnsemblID")]
  tT.EBase = tT.EBase[tT4.EBase, on=c("GeneSymbol", "EnsemblID")]
  
  return(tT.EBase)
}