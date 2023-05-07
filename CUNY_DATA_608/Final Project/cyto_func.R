#' Aggregates Karyotype Abnormalities - Clones
#'
#' Converts feature table into summary table, but splits by clones
#'
#' @param kary.table Features table from `rbindlist.abnormalities`
#'
#' @return kary.table.dcast
#' @examples tbl.dcast = aggregate_tbl.clone(features.tbl)
#' @export
#'


aggregate_tbl.clone <- function(kary.table) {
  if (!is.data.table(kary.table)){
    kary.table.dcast = kary.table
    return(kary.table.dcast)
  }
  if (!is.null(kary.table$subclone)) {
    clone.singular = kary.table[, .N, by = c('clone_no', 'abnormality')][is.na(abnormality) & N == 1, clone_no]
    kary.table.sub = kary.table[is.na(abnormality) & clone_no %in% clone.singular, abnormality := 'Normal']
    kary.table.sub = kary.table[is.na(subclone) & !is.na(abnormality),]
  } else if (!is.null(kary.table$abnormality)) {
    #   kary.table.sub = kary.table[!is.na(abnormality),]
    clone.singular = kary.table[, .N, by = c('clone_no', 'abnormality')][is.na(abnormality) & N == 1, clone_no]
    kary.table.sub = kary.table[is.na(abnormality) & clone_no %in% clone.singular, abnormality := 'Normal']
    kary.table.sub = kary.table[is.na(abnormality), abnormality := 'Normal']
  } else {
    kary.table.sub = kary.table
  }
  if (is.null(kary.table.sub$abnormality)) {
    kary.table.sub[, col_name := "Normal"]
    kary.table.dcast = dcast(kary.table.sub[, list(accession_no, MRN, col_name, clone_no)], MRN + accession_no + clone_no ~ col_name, fun=length)
    if("chromosomal_count" %in% names(kary.table)) {
      kary.table.dcast[, chromosomal_count := toString(unique(kary.table[, (chromosomal_count)]))]
    }
  } else {
    if (is.null(kary.table.sub$count) & !(is.null(kary.table.sub$count1) & is.null(kary.table.sub$count2))) {
      kary.table.sub[, count := paste0(count1, "_", count2)]
    }
    if (c("chrom1") %in% colnames(kary.table.sub)) {
      kary.table.sub[, chrom_var := chrom1]
      kary.table.sub[(abnormality == "marker" | abnormality == "double-minute"), chrom_var := chrom1]
      kary.table.sub[(abnormality == "translocation" & !("chrom3" %in% names(kary.table.sub))), chrom_var := paste0(chrom1, ".", chrom2)]
      kary.table.sub[(abnormality == "translocation" & ("chrom3" %in% names(kary.table.sub))), chrom_var := paste0(chrom1, ".", chrom2, ".", chrom3)]
      kary.table.sub[(abnormality == "Incomplete" | abnormality == "UNKNOWN"), chrom_var := NA]
    } else {
      kary.table.sub[, chrom_var := NA]
      kary.table.sub[abnormality == "marker" | abnormality == "double-minute", chrom_var := count]
      kary.table.sub[abnormality == "translocation" & !("chrom3" %in% names(kary.table.sub)), chrom_var := paste0(chrom1, ".", chrom2)]
      kary.table.sub[abnormality == "translocation" & ("chrom3" %in% names(kary.table.sub)), chrom_var := paste0(chrom1, ".", chrom2, ".", chrom3)]
      kary.table.sub[(abnormality == "Incomplete" | abnormality == "UNKNOWN"), chrom_var := NA]
      
    }
    kary.table.sub[, region_var := ""]
    kary.table.sub[(abnormality == "translocation" | abnormality == "inversion" | abnormality == "duplication")  & length(intersect(names(kary.table.sub), c("region1", "region2"))) == 2, region_var := paste0(region1, ".", region2)]
    # kary.table.sub[(abnormality == "translocation" | abnormality == "inversion" | abnormality == "duplication")  & length(intersect(names(kary.table.sub), c("region1", "region2", "region3"))) == 3,][!is.na(region3), region_var := paste0(region1, ".", region2, ".", region3)]
    if (length(intersect("region3",names(kary.table.sub))) == 1)  {
      # kary.table.sub[(abnormality == "translocation")  & length(intersect(names(kary.table.sub), c("region1", "region2", "region3"))) == 3, region_var := paste0(region1, ".", region2, ".", region3)]
      kary.table.sub[(abnormality == "translocation")  & (!is.na(region3) & region3 != "NA"), region_var := paste0(region1, ".", region2, ".", region3)]
    }
    kary.table.sub[abnormality == "derivative" & length(intersect(names(kary.table.sub), c("region1", "region2"))) == 2, region_var := paste0(region1, ".", region2)]
    kary.table.sub[abnormality == "derivative" & length(intersect(names(kary.table.sub), c("region11", "region12"))) == 2, region_var := paste0(region11, ".", region12)]
    kary.table.sub[abnormality == "Incomplete" | abnormality == "UNKNOWN", region_var := NA]
    kary.table.sub[abnormality == "marker", region_var := ""]
    kary.table.sub[(abnormality %in% c("deletion", "addition", "ring") & length(intersect(names(kary.table.sub), c("region1"))) == 1), region_var := region1]
    kary.table.sub[(length(intersect(names(kary.table.sub), c("region1"))) == 1) & region_var == "", region_var := region1]
    
    kary.table.sub[, col_name := paste(abnormality,chrom_var, region_var, sep = ".")]
    kary.table.dcast = dcast(kary.table.sub[, list(MRN, accession_no, clone_no, col_name)], MRN + accession_no + clone_no ~ col_name, fun=length)
    if("chromosomal_count" %in% names(kary.table)) {
      kary.table.dcast = merge(unique(kary.table[, .(clone_no, chromosomal_count)]), kary.table.dcast, by = 'clone_no')
      # kary.table.dcast[, chromosomal_count := toString(unique(kary.table[, (chromosomal_count)]))]
    }
  }
  if("clone_no" %in% names(kary.table) & is.data.table(kary.table)) {
    kary.table.dcast[, clone.count := length(unique(kary.table[, clone_no]))]
  }
  if("ploidy" %in% names(kary.table) & is.data.table(kary.table)) {
    kary.table.dcast = merge(unique(kary.table[, .(clone_no, ploidy)]), kary.table.dcast, by = 'clone_no')
  }
  if("cp" %in% names(kary.table) & is.data.table(kary.table)) {
    kary.table.dcast = merge(unique(kary.table[, .(clone_no, cp)]), kary.table.dcast, by = 'clone_no')
  }
  if("NA.NA.NA" %in% names(kary.table.dcast) & exists("kary.table.dcast")) {
    kary.table.dcast[,NA.NA.NA := NULL]
    kary.table.dcast[, Normal := T]
  }
  if('Normal.NA.NA' %in% names(kary.table.dcast) & exists("kary.table.dcast")) {
    setnames(kary.table.dcast, 'Normal.NA.NA', 'Normal')
  }
  if('Normal.NA.' %in% names(kary.table.dcast) & exists("kary.table.dcast")) {
    setnames(kary.table.dcast, 'Normal.NA.', 'Normal')
  }
  if("cell_count" %in% names(kary.table) & is.data.table(kary.table)) {
    kary.table.dcast = merge(unique(kary.table[, .(clone_no, cell_count)]), kary.table.dcast, by = 'clone_no')
  }
  return(kary.table.dcast)
}


#' Function to extract cell count of clone
#'
#' @param text.example Modal karyotype
#' @param extract.report.text Extract cell counts from a pathology report text (Default: `TRUE`)
#' @return cell.counts
#' @export
#'
#'

cell_count_total.clone <- function(text.example,  extract.report.text = T)  {
  if(extract.report.text) {
    text.example = as.character(text.example)
    text.example = cleanPath.str (path.str = text.example)
    text.example = text.example[grepl("Modal karyotype[:]|Modal karyotype[(]s[)][:]", text.example)]
    if (length(text.example) == 0) {
      return(NA)
      break
    }
    
    if (grepl("Modal karyotype[:]", text.example)) {
      text.example = gsub(".*Modal karyotype[:]", "", text.example)  #Isolate Karyotype, left
    } else {
      text.example = gsub(".*Modal karyotype[(]s[)][:]", "", text.example)  #Isolate Karyotype, left
    }
    
    text.example = gsub("Cytoscan.*|ANALYSIS.*|FISH.*|nuc.*|TEST.*|DISCLAIMER.*|DIAGNOSTIC.*|e-mailed.*|BONE.*|Note.*|Interim.*|PREPARATION.*|Probe.*|yes.*|ABNORMAL.*|not.*|Int(.)?.*", "", text.example)
    
    text.example = gsub('^/', '', text.example) #remove leading blank
    text.example = gsub('\\/', '', text.example) #remove slash
    text.example = gsub("\\n", "", text.example, perl = T) # Remove Carriage Returns
    text.example = gsub("\\\\n", "", text.example)
    text.example = gsub("idemdel", "idem,del", text.example) # Remove Typos
    text.example = gsub("\\s", "", text.example) # Remove White Spaces
    text.example = gsub("\\?", "", text.example) # Remove Question Mark
    
    
    text.example = unlist(strsplit(text.example, "]")) # Split on square bracket , unlist
    text.example = paste0(text.example, "]") # Add back bracket
    text.example = setdiff(text.example, "]")
    text.example = setdiff(text.example, ".]")
    if (length(text.example) == 0) { # For bad reports :(
      return(NA)
    } else {
      
      lst.text = split(text.example, 1:length(text.example)) # split karyotype into list
      
      cell.counts = gsub("\\[|\\]", "", regmatches(lst.text, gregexpr("\\[.*?\\]", lst.text))) # Extract Cell Counts
      cell.counts = gsub("cp", "", cell.counts) # Swap out CP
      cell.counts = as.integer(cell.counts)
      #  cell.total = sum(cell.counts)
    }
    return(cell.counts)
  } else {
    text.example = gsub('^/', '', text.example) #remove leading blank
    text.example = gsub('\\/', '', text.example) #remove slash
    text.example = gsub("\\n", "", text.example, perl = T) # Remove Carriage Returns
    text.example = gsub("\\\\n", "", text.example)
    text.example = gsub("idemdel", "idem,del", text.example) # Remove Typos
    text.example = gsub("\\s", "", text.example) # Remove White Spaces
    text.example = gsub("\\?", "", text.example) # Remove Question Mark
    
    
    text.example = unlist(strsplit(text.example, "]")) # Split on square bracket , unlist
    text.example = paste0(text.example, "]") # Add back bracket
    text.example = setdiff(text.example, "]")
    text.example = setdiff(text.example, ".]")
    if (length(text.example) == 0) { # For bad reports :(
      return(NA)
    } else {
      
      lst.text = split(text.example, 1:length(text.example)) # split karyotype into list
      
      cell.counts = gsub("\\[|\\]", "", regmatches(lst.text, gregexpr("\\[.*?\\]", lst.text))) # Extract Cell Counts
      cell.counts = gsub("cp", "", cell.counts) # Swap out CP
      cell.counts = as.integer(cell.counts)
      #  cell.total = sum(cell.counts)
    }
    return(cell.counts)
  }
}

#' Function to parse abnormalities per clone
#'
#' @param dt Data Table w cytogenetics
#' @param col.MRN MRN Column
#' @param col.accession.no Accession.No Column
#' @param cell_count_lim Cell count for incomplete flagging
#' @param report.text Report text
#'
#' @return dt.cyto
#' @export
#'
parseKaryotypewrapper.clone = function (dt, col.MRN = "MRN", col.accession.no = "Accession.No",   report.text = "nan", extract.report.text = T) {
  if(extract.report.text == T) {
    dt[, karyotype_included := sapply(get(report.text), karyotype_included)]
    dt.path.cyto.karyotype = dt[karyotype_included == T,]
  } else {
    dt.path.cyto.karyotype = copy(dt)
  }
  lst.tbls = list()
  for (i in 1:nrow(dt.path.cyto.karyotype)) {
    if(extract.report.text == T) {
      text.example = karyotype_parser(dt.path.cyto.karyotype[i, get(report.text)])
    } else {
      text.example = karyotype_parser.non_report(dt.path.cyto.karyotype[i, get(report.text)])
    }
    parsed = flag_function(text.example)
    kary.table = rbindlist.abnormalities(parsed, MRN = dt.path.cyto.karyotype[i,  get(col.MRN)], accession.no = dt.path.cyto.karyotype[i, get(col.accession.no)], cell_count_lim = 0)
    kary.tbl.wide = aggregate_tbl.clone(kary.table)
    #   print(names(kary.tbl.wide))
    #    print(names(kary.table))
    if(length(kary.table) > 1) {
      kary.tbl.wide = merge(kary.tbl.wide, unique(kary.table[,.(clone_no, parent_clone)]), by = c('clone_no')) #attempt
    }
    #    print(paste('Looping', dt.path.cyto.karyotype[i, get(col.accession.no)]))
    lst.tbls[[dt.path.cyto.karyotype[i, get(col.accession.no)]]] = kary.tbl.wide
  }
  # return(lst.tbls)
  lst.tbls = lst.tbls[sapply(lst.tbls, function(x) length(x) > 1)]
  vec.accession.no = names(lst.tbls)
  lst.cell_total = sapply(dt[, get(report.text)], cell_count_total.clone, extract.report.text = extract.report.text)
  names(lst.cell_total) = as.character(dt[, get(col.accession.no)])
  dt.cell.total = melt(lst.cell_total)
  dt.cell.total = as.data.table(dt.cell.total)
  setnames(dt.cell.total, c('cell.total', 'Accession.No'))
  dt.cell.count = dt.cell.total[Accession.No %in% vec.accession.no]
  dt.cell.count[,  clone.ID := 1:.N, by = 'Accession.No']
  
  dt.path.cyto.karyotype[, `:=`(cell_total, sapply(get(report.text), cell_count_total))]
  # dt.cyto = cytogenetic.func(lst.tbls)
  dt.cyto = rbindlist(lst.tbls, fill = T)
  setcolorder(dt.cyto, c("MRN", "accession_no", "chromosomal_count", "clone.count", "ploidy", "cp", 'parent_clone'))
  setnames(dt.cyto, 'accession_no', "Accession.No")
  dt.cyto[is.na(dt.cyto)] <- 0
  dt.cyto[,  clone.ID := 1:.N, by = 'Accession.No']
  cols = c(col.MRN, col.accession.no, "cell_total")
  dt.cyto = merge(dt.cyto, dt.path.cyto.karyotype[, .SD, .SDcols = cols],  by.x = c("MRN", "Accession.No"), by.y = c(col.MRN, col.accession.no),  all.x = T)
  dt.cyto = merge(dt.cyto, dt.cell.count, by = c('Accession.No', 'clone.ID'))
  
  return(dt.cyto)
}

#' @import data.table
NULL

#' Assigns ELN Abnormalities
#'
#' Creates a table with the relevant ELN abnormalities
#'
#' @param dt.cyto A table with parsed abnormalities
#'
#' @return lst.checked_karyotypes.binary
#' @examples dt.path.cytoDMP.NCCN = getCyto.ELN(dt.path.cytoDmp.dx)
#' @export
#'

getCyto.ELN <- function(dt.cyto = dt.cyto) {
  dt.matching.regex = fread(system.file('extdata', 'nlp_parser_fish_matching_2022.txt', package='karyoParser'))
  parser_regex = unlist(dt.matching.regex$parser_regex)
  parser_regex = parser_regex[parser_regex != ""]
  cyto_regex = unlist(dt.matching.regex$leukNLP_Regex)
  
  lst.cyto.abnormalities = list()
  for (i in cyto_regex) {
    j = unlist(dt.matching.regex[leukNLP_Regex == i, "parser_regex"], use.names = F)
    if (j == "") {
      next
    } else {
      x = names(dt.cyto)[grep(j, names(dt.cyto))]
      lst.cyto.abnormalities[[i]] = x
    }
  }
  
  for (i in (1:length(lst.cyto.abnormalities))) {
    abnormality.cols = unlist(lst.cyto.abnormalities[[i]])
    colname = unlist(names(lst.cyto.abnormalities[i]))
    if (length(abnormality.cols) == 0) {
      dt.cyto[, (colname) := FALSE]
    } else {
      dt.cyto[, (colname) :=   rowSums(dt.cyto[, .SD, .SDcols = abnormality.cols])]
    }
  }
  
  # Convert Monosomal to partial
  dt.cyto[del7 == T, `:=` (del7p = T, del7q = T)]
  dt.cyto[del17 == T, `:=` (del17p = T)]
  dt.cyto[del5 == T, `:=` (del5p = T, del5q = T)]
  dt.cyto[del6 == T, `:=` (del6p = T)]
  
  
  
  # Convert to Binary
  dt.cyto = cbind(dt.cyto[, .SD, .SDcols = (setdiff(names(dt.cyto), cyto_regex))],
                  dt.cyto[, .SD, .SDcols = (intersect(names(dt.cyto), cyto_regex))][, lapply(.SD, function(x) x = ifelse(x == 0 | x == F, F, T))])
  dt.cyto[, cyto_label :=  as.character(cytogenetics)]
  dt.cyto[del5q == T, cyto_label := "del5q"]
  dt.cyto[del17 == T, cyto_label := "del17"]
  dt.cyto[del7 == T, cyto_label := "del7"]
  dt.cyto[abn17p == T, cyto_label := "abn(17p)"]
  dt.cyto[tv.3 == T, cyto_label := 't(3q26.2;v)']
  dt.cyto[tv.11 == T, cyto_label := "t(v;11q23)"]
  dt.cyto[t6.9 == T, cyto_label := "t6.9"]
  dt.cyto[inv3 == T, cyto_label := "inv(3)"]
  dt.cyto[cytogenetics == "Monosomal", cyto_label :=  as.character(cytogenetics)]
  dt.cyto[cytogenetics == "Complex", cyto_label :=  as.character(cytogenetics)]
  dt.cyto[cytogenetics %in% c("CBF", "APL"), cyto_label :=  as.character(cytogenetics)]
  dt.cyto[cytogenetics == 'Unknown', cyto_label := 'Unknown']
  
  dt.cyto[, cyto_label := as.character(cyto_label)]
  return(dt.cyto)
}


#' Karyotype Parser Wrapper Function
#'
#' Creates a table with all abnormalities, cytogenetics category, and cell count
#'
#' @param dt A table with report text, MRN, and accession number
#' @param col.MRN Column containing MRN
#' @param col.accession.no Column containing Accession No
#' @param cell_count_lim Limit for cell count to be parsed and incorporataed into karyotype (defaults at 1)
#' @param report.text Full report text column (default is nan)
#' @param AML Is the parsing AML specific? Defaulted to `TRUE`
#' @param congenital.abnormality Do congential abnormalities count towards totals? Default: `FALSE`
#' @param ELN `2017` or `2022` (Default: `2022`). Affects Complex Karyotype definition
#'
#' @return dt.cyto
#' @examples dt.cyto = parseKaryotypewrapper(dt)
#' @export
#'
#'

parseKaryotypewrapper <- function(dt, col.MRN = "MRN", col.accession.no = "Accession.No", cell_count_lim = 1, report.text = "nan", AML = T, extract.report.text = T, congenital.abnormality = F, ELN = '2022') {
  
  if(extract.report.text == T) {
    dt[, karyotype_included := sapply(get(report.text), karyotype_included)]
    dt.path.cyto.karyotype = dt[karyotype_included == T,]
  } else {
    dt.path.cyto.karyotype = copy(dt)
  }
  
  
  lst.tbls = list()
  for (i in 1:nrow(dt.path.cyto.karyotype)) {
    if(extract.report.text == T) {
      text.example = karyotype_parser(dt.path.cyto.karyotype[i, get(report.text)])
    } else {
      text.example = karyotype_parser.non_report(dt.path.cyto.karyotype[i, get(report.text)])
    }
    parsed = flag_function(text.example)
    kary.table = rbindlist.abnormalities(parsed, MRN = dt.path.cyto.karyotype[i, get(col.MRN)], accession.no = dt.path.cyto.karyotype[i, get(col.accession.no)], cell_count_lim = (cell_count_lim))
    kary.tbl.wide = aggregate_tbl(kary.table)
    lst.tbls[[i]] = kary.tbl.wide
  }
  
  
  
  ## Remove Karyotypes with a cell count of 1
  lst.tbls = lst.tbls[sapply(lst.tbls, function(x) length(x) > 1)]
  
  ## Summing Cell Counts Per Karyotype
  dt.path.cyto.karyotype[, cell_total := sapply(get(report.text), cell_count_total, extract.report.text = extract.report.text)]
  
  ## Determine Cytogenetics
  if (AML == T ){
    dt.cyto = cytogenetic.func(lst.tbls, congenital.abnormality = congenital.abnormality, ELN = ELN)
  } else {
    dt.cyto = cytogenetic.func.non.AML(lst.tbls, congenital.abnormality = congenital.abnormality)
  }
  cols = c(col.MRN, col.accession.no, "cell_total")
  dt.cyto = merge(dt.cyto, dt.path.cyto.karyotype[, .SD, .SDcols = cols], by.x = c("MRN", "Accession.No"), by.y = c(col.MRN, col.accession.no), all.x = T)
  
  return(dt.cyto)
}

#' @import data.table
NULL

#' Create Cytogenetics Abnormality Table
#'
#' Produces a table of abnormalities for the AML Risk assignment function
#'
#' @param dt.parsed A table with parsed abnormalities
#' @param ELN Risk year (takes either `2017` or `2022`)
#'
#' @return dt.risk.cyto
#' @examples dt.risk.cyto = createAbnormality.tbl(dt.parsed)
#'
#'


createAbnormality.tbl <-  function(dt.parsed, ELN = '2017') {
  if (ELN == '2017') {
    lbl.risk.cyto.revised = fread(system.file('extdata', 'risk_cyto_ELN_2017.txt', package='karyoParser'))
  } else if (ELN == '2022') {
    lbl.risk.cyto.revised = fread(system.file('extdata', 'risk_cyto_ELN_2022.txt', package='karyoParser'))
  } else {
    return('Please specify ELN year')
  }
  
  abnormality.tbl = lbl.risk.cyto.revised[feature != "cytogenetics",]
  feature = c()
  for (i in 1:nrow(abnormality.tbl)){
    cyto.risk.cols = colnames(dt.parsed)[grepl(abnormality.tbl[i,feature], colnames(dt.parsed))]
    names(cyto.risk.cols) = rep(abnormality.tbl$risk[i], length(cyto.risk.cols))
    feature = c(feature, cyto.risk.cols)
  }
  
  dt.risk.cyto = as.data.table(feature, keep.rownames = "risk")
  dt.risk.cyto[, value := T]
  dt.risk.cyto = rbind(lbl.risk.cyto.revised[feature == "cytogenetics",], dt.risk.cyto)
  dt.risk.cyto = unique(dt.risk.cyto)
  
  cyto.risk.cols = c("cytogenetics", feature)
  names(cyto.risk.cols) = NULL
  return(dt.risk.cyto)
}


#' Create Cytogenetics Abnormality Column Vec
#'
#' Produces a vector of abnormalities for the AML Risk assignment function
#'
#' @param dt.parsed A table with parsed abnormalities
#' @param ELN Risk year (takes either `2017` or `2022`)

#' @return cyto.risk.cols
#' @examples createAbnormality.vec(dt.parsed)
#'
#'


createAbnormality.vec <-  function(dt.parsed, ELN = '2017') {
  if (ELN == '2017') {
    lbl.risk.cyto.revised = fread(system.file('extdata', 'risk_cyto_ELN_2017.txt', package='karyoParser'))
  } else if (ELN == '2022') {
    lbl.risk.cyto.revised = fread(system.file('extdata', 'risk_cyto_ELN_2022.txt', package='karyoParser'))
  } else {
    return('Please specify ELN year')
  }
  abnormality.tbl = lbl.risk.cyto.revised[feature != "cytogenetics",]
  feature = c()
  for (i in 1:nrow(abnormality.tbl)){
    cyto.risk.cols = colnames(dt.parsed)[grepl(abnormality.tbl[i,feature], colnames(dt.parsed))]
    names(cyto.risk.cols) = rep(abnormality.tbl$risk[i], length(cyto.risk.cols))
    feature = c(feature, cyto.risk.cols)
  }
  
  dt.risk.cyto = as.data.table(feature, keep.rownames = "risk")
  dt.risk.cyto[, value := T]
  dt.risk.cyto = rbind(lbl.risk.cyto.revised[feature == "cytogenetics",], dt.risk.cyto)
  dt.risk.cyto = unique(dt.risk.cyto)
  
  cyto.risk.cols = c("cytogenetics", feature)
  names(cyto.risk.cols) = NULL
  return(cyto.risk.cols)
}


#' Produce non Intermediate Abnormalities (Vec)
#'
#'
#' @param dt.parsed Produces a vector of abnormalities that can determine cytogenetics w/o dmp
#' @param ELN Risk year (takes either `2017` or `2022`)
#'
#' @return cyto.risk.cols
#' @examples cytoAbnormalities.vec(dt.parsed)
#'
#'


cytoAbnormalities.vec.nonIntermediate <-  function(dt.parsed, ELN = '2017') {
  if (ELN == '2017') {
    lbl.risk.cyto.revised = fread(system.file('extdata', 'risk_cyto_ELN_2017.txt', package='karyoParser'))
  } else if (ELN == '2022') {
    lbl.risk.cyto.revised = fread(system.file('extdata', 'risk_cyto_ELN_2022.txt', package='karyoParser'))
  } else {
    return('Please specify ELN year')
  }
  
  abnormality.tbl = lbl.risk.cyto.revised[feature != "cytogenetics" & risk != "Intermediate",]
  feature = c()
  for (i in 1:nrow(abnormality.tbl)){
    cyto.risk.cols = colnames(dt.parsed)[grepl(abnormality.tbl[i,feature], colnames(dt.parsed))]
    names(cyto.risk.cols) = rep(abnormality.tbl$risk[i], length(cyto.risk.cols))
    feature = c(feature, cyto.risk.cols)
  }
  
  dt.risk.cyto = as.data.table(feature, keep.rownames = "risk")
  dt.risk.cyto[, value := T]
  dt.risk.cyto = rbind(lbl.risk.cyto.revised[feature == "cytogenetics",], dt.risk.cyto)
  dt.risk.cyto = unique(dt.risk.cyto)
  
  cyto.risk.cols = c(feature)
  names(cyto.risk.cols) = NULL
  return(cyto.risk.cols)
}

#' @import data.table
NULL

#' Assigns Cytogenetics Categories
#'
#' Assigns feature type to each karyotype abnormality
#'
#' @param lst.kary.wide A list of tables from the `aggregate.tbl` function.
#' @param congenital.abnormality Do congenital abnormalities count toward abnormality counts? Defaulted to `FALSE`
#' @param ELN `2017` or `2022` (Default: `2022`). Affects Complex Karyotype definition
#'
#' @return lst.checked_karyotypes.binary
#' @examples lst.checked_karyotypes.binary = cytogenetic.func(lst.kary.wide)
#' @export
#'


cytogenetic.func <- function(lst.kary.wide, congenital.abnormality = F, ELN = '2022') {
  
  lst.checked_karyotypes.rbind = rbindlist(lst.kary.wide, fill = T)
  setcolorder(lst.checked_karyotypes.rbind, c("MRN", "accession_no", "chromosomal_count", "clone.count", "ploidy", "cp"))
  
  lst.checked_karyotypes.binary = as.data.table(apply(lst.checked_karyotypes.rbind[, 7:ncol(lst.checked_karyotypes.rbind)], 2, function(x) x <- ifelse(is.na(x), F, T)))
  
  lst.checked_karyotypes.binary[, MRN := lst.checked_karyotypes.rbind$MRN][, accession_no := lst.checked_karyotypes.rbind$accession_no][, chromosomal_count := lst.checked_karyotypes.rbind$chromosomal_count][, clone.count := lst.checked_karyotypes.rbind$clone.count][, ploidy := lst.checked_karyotypes.rbind$ploidy][, cp := lst.checked_karyotypes.rbind$cp]
  
  setcolorder(lst.checked_karyotypes.binary, c("MRN", "accession_no", "chromosomal_count", "clone.count", "ploidy", "cp"))
  
  
  full.deletion.cols = names(lst.checked_karyotypes.binary)[setdiff(intersect(grep("FULL", names(lst.checked_karyotypes.binary)),
                                                                              grep("deletion", names(lst.checked_karyotypes.binary))),
                                                                    grep("der", names(lst.checked_karyotypes.binary)))]
  
  full.deletion.cols = setdiff(full.deletion.cols,  c("deletion.X.FULL", "deletion.Y.FULL")) # Remove X and Y
  
  
  # Tetraploid Evaluation
  lst.checked_karyotypes.binary[, polyploid := ifelse(grepl("^8[0-9]|^9[0-9]|88|92", chromosomal_count, perl = T), T, F)]
  lst.checked_karyotypes.binary[polyploid == F & grepl("3n|4n", ploidy), polyploid := T]
  
  # Marker Chromosomes
  vec.marker.cols = names(lst.checked_karyotypes.binary)[grep("marker", names(lst.checked_karyotypes.binary))]
  
  # Congenital Abnormalities
  vec.congential.cols = names(lst.checked_karyotypes.binary)[grep("congenital", names(lst.checked_karyotypes.binary))]
  
  # List full addition cols
  full.addition.cols = names(lst.checked_karyotypes.binary)[grep("addition.(.)+.FULL", names(lst.checked_karyotypes.binary))]
  
  # List ring cols
  ring.cols = names(lst.checked_karyotypes.binary)[grep("ring", names(lst.checked_karyotypes.binary))]
  
  # List double-minute cols
  dm.cols = names(lst.checked_karyotypes.binary)[grep("double-minute", names(lst.checked_karyotypes.binary))]
  
  # Check for Full Deletions
  lst.checked_karyotypes.binary$full_deletion <- unlist(apply(lst.checked_karyotypes.binary[, ..full.deletion.cols], 1, function(x) x <- ifelse(any(x), 1, NA)))
  lst.checked_karyotypes.binary$full_deletion  = ifelse(is.na(lst.checked_karyotypes.binary$full_deletion), 0, lst.checked_karyotypes.binary$full_deletion)
  
  # Return count of full deletions
  lst.checked_karyotypes.binary$total_full_deletions = rowSums(lst.checked_karyotypes.binary[, ..full.deletion.cols])
  
  # Return count of full additions
  lst.checked_karyotypes.binary$total_full_additions = rowSums(lst.checked_karyotypes.binary[, ..full.addition.cols])
  
  # CBF Features
  cbf.vec = c(names(lst.checked_karyotypes.binary)[grep("inversion.16", names(lst.checked_karyotypes.binary))],
              names(lst.checked_karyotypes.binary)[grep("translocation.8.21|translocation.([0-9]+[X|Y])\\.8.21", names(lst.checked_karyotypes.binary))],
              names(lst.checked_karyotypes.binary)[grep("derivative.[0-9]+.8.21", names(lst.checked_karyotypes.binary))]
  )
  
  # APL Features
  apl.vec = c(names(lst.checked_karyotypes.binary)[grep("translocation.15.17|translocation.([0-9]+[X|Y])\\.15.17", names(lst.checked_karyotypes.binary))],
              names(lst.checked_karyotypes.binary)[grep("derivative.[0-9]+.15.17", names(lst.checked_karyotypes.binary))])
  
  # Abnormality Columns
  if (congenital.abnormality == F) {
    abnormalities.cols = setdiff(names(lst.checked_karyotypes.binary), c("MRN", "accession_no", "chromosomal_count", "clone.count", "ploidy", "cp", "full_deletion", "deletion.X.FULL", "deletion.Y.FULL", "MK.parsed", "MK", "CBF", "CK (3 or more)", "Incomplete.NA.NA", cbf.vec, apl.vec, vec.congential.cols, "complex", "Normal", "complex_BINARY", "total_full_deletions", "total_full_additions", full.deletion.cols, names(lst.checked_karyotypes.binary)[grep("marker", names(lst.checked_karyotypes.binary))])) # WHO NOT INCLUDED IN MONOSOMY
  } else {
    abnormalities.cols = setdiff(names(lst.checked_karyotypes.binary), c("MRN", "accession_no", "chromosomal_count", "clone.count", "ploidy", "cp", "full_deletion", "deletion.X.FULL", "deletion.Y.FULL", "MK.parsed", "MK", "CBF", "CK (3 or more)", "Incomplete.NA.NA", cbf.vec, apl.vec, "complex", "Normal", "complex_BINARY", "total_full_deletions", "total_full_additions", full.deletion.cols, names(lst.checked_karyotypes.binary)[grep("marker", names(lst.checked_karyotypes.binary))])) # WHO NOT INCLUDED IN MONOSOMY   # Remove marker from MK/Complex abnormalities
    
  }
  
  monosomy.abnormality.cols = setdiff(abnormalities.cols, vec.marker.cols) # Remove marker from monosomy abnormalities
  monosomy.abnormality.cols = setdiff(monosomy.abnormality.cols, dm.cols)# Remove double minute from monosomy abnormalities
  monosomy.abnormality.cols = setdiff(monosomy.abnormality.cols, full.addition.cols) # Remove full addition from monosomy
  monosomy.abnormality.cols = setdiff(monosomy.abnormality.cols, ring.cols)   # Remove ring from monosomy abnormalities
  monosomy.abnormality.cols = setdiff(monosomy.abnormality.cols, "ploidy") # Remove ploidy from monosomal abnormalities
  monosomy.abnormality.cols = setdiff(monosomy.abnormality.cols,  c("deletion.X.FULL", "deletion.Y.FULL")) # Remove X and Y
  
  lst.checked_karyotypes.binary[full_deletion == 1, MK_additional_abnormalities:= unlist(apply(lst.checked_karyotypes.binary[full_deletion== 1, .SD, .SDcols = monosomy.abnormality.cols] , 1, function(x) x <- sum(x, na.rm = T)))]
  lst.checked_karyotypes.binary[, MK_additional_abnormalities:= ifelse(is.na(MK_additional_abnormalities), 0, MK_additional_abnormalities)]
  
  # t(9;11) Columns
  t9.11.cols = names(lst.checked_karyotypes.binary)[grep("translocation.9.11", names(lst.checked_karyotypes.binary))]
  
  # Monosomy evaluation
  #  lst.checked_karyotypes.binary[, MK_count := full_deletion + MK_additional_abnormalities]
  lst.checked_karyotypes.binary[, MK_count := total_full_deletions + MK_additional_abnormalities]
  lst.checked_karyotypes.binary[, MK := ifelse(MK_count > 1, 1, 0)]
  
  who.vec = c(names(lst.checked_karyotypes.binary)[grep("translocation\\.9\\.11|derivative\\.[0-9]+\\.9\\.11", names(lst.checked_karyotypes.binary))],
              names(lst.checked_karyotypes.binary)[grep("translocation.6.9|translocation.([0-9]+[X|Y])\\.6.9|derivative.[0-9]+.6.9", names(lst.checked_karyotypes.binary))],
              names(lst.checked_karyotypes.binary)[grep("inversion.3", names(lst.checked_karyotypes.binary))],
              names(lst.checked_karyotypes.binary)[grep("inversion.16", names(lst.checked_karyotypes.binary))],
              names(lst.checked_karyotypes.binary)[grep("translocation.[0-9]+\\.11\\.[p|q]\\d+\\.[^p]", names(lst.checked_karyotypes.binary))],
              names(lst.checked_karyotypes.binary)[grep("translocation.11\\.[0-9]+\\.[^p]|translocation.([0-9]+[X|Y])\\.11", names(lst.checked_karyotypes.binary))],
              names(lst.checked_karyotypes.binary)[grep("translocation.3.3", names(lst.checked_karyotypes.binary))],
              names(lst.checked_karyotypes.binary)[grep("translocation.16.16", names(lst.checked_karyotypes.binary))])
  
  # Complex Evaluation --> WHO counts for abnormlaities
  if (congenital.abnormality == F) {
    lst.checked_karyotypes.binary$abnormality_total <- unlist(apply(lst.checked_karyotypes.binary[, .SD, .SDcols = setdiff(colnames(lst.checked_karyotypes.binary), c("MRN", "Normal", "accession_no",  "clone.count", "chromosomal_count", "ploidy", "cp", "full_deletion", "total_full_deletions", "MK_count", "MK", "CBF", "APL", "ck_3_or_more", "Incomplete.NA.NA", "total_full_additions", cbf.vec, apl.vec, who.vec, vec.congential.cols))], 1, function(x) x <- sum(x, na.rm = T)))
  } else {
    lst.checked_karyotypes.binary$abnormality_total <- unlist(apply(lst.checked_karyotypes.binary[, .SD, .SDcols = setdiff(colnames(lst.checked_karyotypes.binary), c("MRN", "Normal", "accession_no", "clone.count", "chromosomal_count",  "ploidy", "cp", "full_deletion", "total_full_deletions", "MK_count", "MK", "CBF", "APL", "ck_3_or_more", "Incomplete.NA.NA", "total_full_additions", cbf.vec, apl.vec, who.vec))], 1, function(x) x <- sum(x, na.rm = T)))
  }
  
  
  # Add a Raw Abnormality Total (includes WHO, etc)
  if (congenital.abnormality == F) {
    lst.checked_karyotypes.binary$raw_abnormality_total <- unlist(apply(lst.checked_karyotypes.binary[, .SD, .SDcols = setdiff(colnames(lst.checked_karyotypes.binary), c("MRN", "Normal", "accession_no",  "clone.count", "chromosomal_count", "ploidy", "cp", "full_deletion", "total_full_deletions", "MK_count", "MK", "CBF", "APL", "ck_3_or_more", "Incomplete.NA.NA", "total_full_additions", "abnormality_total", "MK_additional_abnormalities", vec.congential.cols))], 1, function(x) x <- sum(x, na.rm = T)))
  } else {
    lst.checked_karyotypes.binary$raw_abnormality_total <- unlist(apply(lst.checked_karyotypes.binary[, .SD, .SDcols = setdiff(colnames(lst.checked_karyotypes.binary), c("MRN", "Normal", "accession_no", "clone.count", "chromosomal_count",  "ploidy", "cp", "full_deletion", "total_full_deletions", "MK_count", "MK", "CBF", "APL", "ck_3_or_more", "Incomplete.NA.NA", "total_full_additions",  "abnormality_total", "MK_additional_abnormalities"))], 1, function(x) x <- sum(x, na.rm = T)))
  }
  
  lst.checked_karyotypes.binary$complex_BINARY  = ifelse(lst.checked_karyotypes.binary$abnormality_total > 2, 1, 0)
  if(ELN == '2022') {
    lst.checked_karyotypes.binary$complex_BINARY  = ifelse(lst.checked_karyotypes.binary$abnormality_total > 2 & lst.checked_karyotypes.binary$abnormality_total != lst.checked_karyotypes.binary$total_full_additions, 1, 0)
  }
  lst.checked_karyotypes.binary$CBF <- unlist(apply(lst.checked_karyotypes.binary[, ..cbf.vec], 1, function(x) x <- ifelse(any(x), 1, 0)))
  lst.checked_karyotypes.binary$APL <- unlist(apply(lst.checked_karyotypes.binary[, ..apl.vec], 1, function(x) x <- ifelse(any(x), 1, 0)))
  lst.checked_karyotypes.binary$t9.11<- unlist(apply(lst.checked_karyotypes.binary[, ..t9.11.cols], 1, function(x) x <- ifelse(any(x), 1, 0)))
  lst.checked_karyotypes.binary$congenital <- unlist(apply(lst.checked_karyotypes.binary[, ..vec.congential.cols], 1, function(x) x <- ifelse(any(x), 1, 0)))
  
  
  # Accounting for NAs
  lst.checked_karyotypes.binary$MK = ifelse(is.na(lst.checked_karyotypes.binary$MK), 0, lst.checked_karyotypes.binary$MK)
  lst.checked_karyotypes.binary$complex_BINARY = ifelse(is.na(lst.checked_karyotypes.binary$complex_BINARY), 0, lst.checked_karyotypes.binary$complex_BINARY)
  lst.checked_karyotypes.binary$CBF = ifelse(is.na(lst.checked_karyotypes.binary$CBF), 0 , lst.checked_karyotypes.binary$CBF)
  lst.checked_karyotypes.binary$APL  = ifelse(is.na( lst.checked_karyotypes.binary$APL), 0,  lst.checked_karyotypes.binary$APL)
  lst.checked_karyotypes.binary$t9.11 = ifelse(is.na( lst.checked_karyotypes.binary$t9.11), 0,  lst.checked_karyotypes.binary$t9.11)
  lst.checked_karyotypes.binary$congenital = ifelse(is.na( lst.checked_karyotypes.binary$congenital), 0,  lst.checked_karyotypes.binary$congenital)
  
  # Labels
  lst.checked_karyotypes.binary[, cytogenetics := 'OND']
  if ("Normal" %in% names(lst.checked_karyotypes.binary)) {
    lst.checked_karyotypes.binary[Normal == T & abnormality_total == 0, cytogenetics := 'Normal']
  }
  if(congenital.abnormality == F) {
    lst.checked_karyotypes.binary[congenital == 1 & abnormality_total == 0, cytogenetics := 'Normal']
  }
  lst.checked_karyotypes.binary[t9.11 == 1, cytogenetics := 't9.11']
  lst.checked_karyotypes.binary[complex_BINARY == 1, cytogenetics := 'Complex']
  lst.checked_karyotypes.binary[MK == 1, cytogenetics := 'Monosomal']
  lst.checked_karyotypes.binary[CBF == 1, cytogenetics := 'CBF']
  lst.checked_karyotypes.binary[APL == 1, cytogenetics := 'APL']
  
  lst.checked_karyotypes.binary$cytogenetics = factor( lst.checked_karyotypes.binary$cytogenetics , levels = c("CBF", "APL", "Monosomal", "Complex", "t9.11", "Normal", "OND"), labels = c("CBF", "APL", "Monosomal", "Complex", "t9.11", "Normal", "OND"))
  #  lst.checked_karyotypes.binary[!(cytogenetics %in% c("CBF", "APL")) & cell_total < 20, cytogenetics := factor("Unknown")]
  
  setnames(lst.checked_karyotypes.binary, "accession_no", "Accession.No")
  setcolorder(lst.checked_karyotypes.binary, c("MRN", "Accession.No", "chromosomal_count",  "ploidy", "cp", "clone.count"))
  
  return(lst.checked_karyotypes.binary)
  
}


#' Creates Wide Table with Tokenized Cytogenetic Abnormalities
#'
#' Assigns feature type to each karyotype abnormality
#'
#' @param lst.kary.wide A list of tables from the `aggregate.tbl` function.
#' @param congenital.abnormality Do congenital abnormalities count as abnormalities? Defaulted to `FALSE`
#' @return lst.checked_karyotypes.binary
#' @examples lst.checked_karyotypes.binary = cytogenetic.func(lst.kary.wide)
#' @export
#'


cytogenetic.func.non.AML <- function(lst.kary.wide, congenital.abnormality = F) {
  
  lst.checked_karyotypes.rbind = rbindlist(lst.kary.wide, fill = T)
  setcolorder(lst.checked_karyotypes.rbind, c("MRN", "accession_no", "chromosomal_count", "clone.count", "ploidy", "cp"))
  
  lst.checked_karyotypes.binary = as.data.table(apply(lst.checked_karyotypes.rbind[, 7:ncol(lst.checked_karyotypes.rbind)], 2, function(x) x <- ifelse(is.na(x), F, T)))
  
  lst.checked_karyotypes.binary[, MRN := lst.checked_karyotypes.rbind$MRN][, accession_no := lst.checked_karyotypes.rbind$accession_no][, chromosomal_count := lst.checked_karyotypes.rbind$chromosomal_count][, clone.count := lst.checked_karyotypes.rbind$clone.count][, ploidy := lst.checked_karyotypes.rbind$ploidy]
  
  setcolorder(lst.checked_karyotypes.binary, c("MRN", "accession_no", "chromosomal_count", "clone.count", "ploidy", "cp"))
  
  
  #  full.deletion.cols = names(lst.checked_karyotypes.binary)[setdiff(intersect(grep("FULL", names(lst.checked_karyotypes.binary)),
  #                                                                              grep("deletion", names(lst.checked_karyotypes.binary))),
  #                                                                    grep("der", names(lst.checked_karyotypes.binary)))]
  
  #  full.deletion.cols = setdiff(full.deletion.cols,  c("deletion.X.FULL", "deletion.Y.FULL")) # Remove X and Y
  
  # Tetraploid Evaluation
  lst.checked_karyotypes.binary[, polyploid := ifelse(grepl("^8[0-9]|^9[0-9]|88|92", chromosomal_count, perl = T), T, F)]
  lst.checked_karyotypes.binary[polyploid == F & grepl("3n|4n", ploidy), polyploid := T]
  
  vec.congential.cols = names(lst.checked_karyotypes.binary)[grep("congenital", names(lst.checked_karyotypes.binary))]
  
  # Add a Raw Abnormality Total (includes WHO, etc)
  if (congenital.abnormality == F) {
    lst.checked_karyotypes.binary$raw_abnormality_total <- unlist(apply(lst.checked_karyotypes.binary[, .SD, .SDcols = setdiff(colnames(lst.checked_karyotypes.binary), c("MRN", "Normal", "accession_no",  "clone.count", "chromosomal_count", "ploidy", "cp", "full_deletion", "total_full_deletions", "MK_count", "MK", "CBF", "APL", "ck_3_or_more", "Incomplete.NA.NA", "abnormality_total", "MK_additional_abnormalities",  vec.congential.cols))], 1, function(x) x <- sum(x, na.rm = T)))
  } else {
    lst.checked_karyotypes.binary$raw_abnormality_total <- unlist(apply(lst.checked_karyotypes.binary[, .SD, .SDcols = setdiff(colnames(lst.checked_karyotypes.binary), c("MRN", "Normal", "accession_no", "clone.count", "chromosomal_count",  "ploidy","cp", "full_deletion", "total_full_deletions", "abnormality_total", "MK_count", "MK_additional_abnormalities", "MK", "CBF", "APL", "ck_3_or_more", "Incomplete.NA.NA"))], 1, function(x) x <- sum(x, na.rm = T)))
  }
  
  # Marker Chromosomes
  #  vec.marker.cols = names(lst.checked_karyotypes.binary)[grep("marker", names(lst.checked_karyotypes.binary))]
  
  # Congenital Abnormalities
  #  vec.congential.cols = names(lst.checked_karyotypes.binary)[grep("congenital", names(lst.checked_karyotypes.binary))]
  
  # List full addition cols
  #  full.addition.cols = names(lst.checked_karyotypes.binary)[grep("addition.(.)+.FULL", names(lst.checked_karyotypes.binary))]
  
  # List ring cols
  #  ring.cols = names(lst.checked_karyotypes.binary)[grep("ring", names(lst.checked_karyotypes.binary))]
  
  # List double-minute cols
  #  dm.cols = names(lst.checked_karyotypes.binary)[grep("double-minute", names(lst.checked_karyotypes.binary))]
  
  # Check for Full Deletions
  #  lst.checked_karyotypes.binary$full_deletion <- unlist(apply(lst.checked_karyotypes.binary[, ..full.deletion.cols], 1, function(x) x <- ifelse(any(x), 1, NA)))
  #  lst.checked_karyotypes.binary$full_deletion  = ifelse(is.na(lst.checked_karyotypes.binary$full_deletion), 0, lst.checked_karyotypes.binary$full_deletion)
  
  # Return count of full deletions
  #  lst.checked_karyotypes.binary$total_full_deletions = rowSums(lst.checked_karyotypes.binary[, ..full.deletion.cols])
  
  # CBF Features
  #  cbf.vec = c(names(lst.checked_karyotypes.binary)[grep("inversion.16", names(lst.checked_karyotypes.binary))],
  #              names(lst.checked_karyotypes.binary)[grep("translocation.8.21|translocation.([0-9]+[X|Y])\\.8.21", names(lst.checked_karyotypes.binary))],
  #              names(lst.checked_karyotypes.binary)[grep("derivative.[0-9]+.8.21", names(lst.checked_karyotypes.binary))]
  #  )
  
  # APL Features
  #  apl.vec = c(names(lst.checked_karyotypes.binary)[grep("translocation.15.17|translocation.([0-9]+[X|Y])\\.15.17", names(lst.checked_karyotypes.binary))],
  #              names(lst.checked_karyotypes.binary)[grep("derivative.[0-9]+.15.17", names(lst.checked_karyotypes.binary))])
  
  # Abnormality Columns
  #  if (congenital.abnormality == F) {
  #    abnormalities.cols = setdiff(names(lst.checked_karyotypes.binary), c("MRN", "accession_no", "chromosomal_count", "clone.count", "ploidy","full_deletion", "deletion.X.FULL", "deletion.Y.FULL", "MK.parsed", "MK", "CBF", "CK (3 or more)", "Incomplete.NA.NA", cbf.vec, apl.vec, vec.congential.cols, "complex", "Normal", "complex_BINARY", "total_full_deletions", full.deletion.cols, names(lst.checked_karyotypes.binary)[grep("marker", names(lst.checked_karyotypes.binary))])) # WHO NOT INCLUDED IN MONOSOMY
  #  } else {
  #    abnormalities.cols = setdiff(names(lst.checked_karyotypes.binary), c("MRN", "accession_no", "chromosomal_count", "clone.count", "ploidy", "full_deletion", "deletion.X.FULL", "deletion.Y.FULL", "MK.parsed", "MK", "CBF", "CK (3 or more)", "Incomplete.NA.NA", cbf.vec, apl.vec, "complex", "Normal", "complex_BINARY", "total_full_deletions", full.deletion.cols, names(lst.checked_karyotypes.binary)[grep("marker", names(lst.checked_karyotypes.binary))])) # WHO NOT INCLUDED IN MONOSOMY   # Remove marker from MK/Complex abnormalities
  
  #  }
  
  #  monosomy.abnormality.cols = setdiff(abnormalities.cols, vec.marker.cols) # Remove marker from monosomy abnormalities
  #  monosomy.abnormality.cols = setdiff(monosomy.abnormality.cols, dm.cols)# Remove double minute from monosomy abnormalities
  #  monosomy.abnormality.cols = setdiff(monosomy.abnormality.cols, full.addition.cols) # Remove full addition from monosomy
  #  monosomy.abnormality.cols = setdiff(monosomy.abnormality.cols, ring.cols)   # Remove ring from monosomy abnormalities
  #  monosomy.abnormality.cols = setdiff(monosomy.abnormality.cols, "ploidy") # Remove ploidy from monosomal abnormalities
  #  monosomy.abnormality.cols = setdiff(monosomy.abnormality.cols,  c("deletion.X.FULL", "deletion.Y.FULL")) # Remove X and Y
  
  #  lst.checked_karyotypes.binary[full_deletion == 1, MK_additional_abnormalities:= unlist(apply(lst.checked_karyotypes.binary[full_deletion== 1, .SD, .SDcols = monosomy.abnormality.cols] , 1, function(x) x <- sum(x, na.rm = T)))]
  #  lst.checked_karyotypes.binary[, MK_additional_abnormalities:= ifelse(is.na(MK_additional_abnormalities), 0, MK_additional_abnormalities)]
  
  # t(9;11) Columns
  #  t9.11.cols = names(lst.checked_karyotypes.binary)[grep("translocation.9.11", names(lst.checked_karyotypes.binary))]
  
  # Monosomy evaluation
  #  lst.checked_karyotypes.binary[, MK_count := full_deletion + MK_additional_abnormalities]
  #lst.checked_karyotypes.binary[, MK_count := total_full_deletions + MK_additional_abnormalities]
  #  lst.checked_karyotypes.binary[, MK := ifelse(MK_count > 1, 1, 0)]
  
  #  who.vec = c(names(lst.checked_karyotypes.binary)[grep("translocation.9.11", names(lst.checked_karyotypes.binary))],
  #              names(lst.checked_karyotypes.binary)[grep("translocation.6.9|translocation.([0-9]+[X|Y])\\.6.9", names(lst.checked_karyotypes.binary))],
  #              names(lst.checked_karyotypes.binary)[grep("inversion.3", names(lst.checked_karyotypes.binary))],
  #              names(lst.checked_karyotypes.binary)[grep("inversion.16", names(lst.checked_karyotypes.binary))],
  #              names(lst.checked_karyotypes.binary)[grep("translocation.[0-9]+\\.11\\.[p|q]\\d+\\.[^p]", names(lst.checked_karyotypes.binary))],
  #              names(lst.checked_karyotypes.binary)[grep("translocation.11\\.[0-9]+\\.[^p]|translocation.([0-9]+[X|Y])\\.11", names(lst.checked_karyotypes.binary))],
  #              names(lst.checked_karyotypes.binary)[grep("translocation.3.3", names(lst.checked_karyotypes.binary))],
  #              names(lst.checked_karyotypes.binary)[grep("translocation.16.16", names(lst.checked_karyotypes.binary))])
  
  # Complex Evaluation --> WHO counts
  #  if (congenital.abnormality == F) {
  #    lst.checked_karyotypes.binary$abnormality_total <- unlist(apply(lst.checked_karyotypes.binary[, .SD, .SDcols = setdiff(colnames(lst.checked_karyotypes.binary), c("MRN", "Normal", "accession_no",  "clone.count", "chromosomal_count", "ploidy", "full_deletion", "total_full_deletions", "MK_count", "MK", "CBF", "APL", "ck_3_or_more", "Incomplete.NA.NA", cbf.vec, apl.vec, who.vec, vec.congential.cols))], 1, function(x) x <- sum(x, na.rm = T)))
  #  } else {
  #    lst.checked_karyotypes.binary$abnormality_total <- unlist(apply(lst.checked_karyotypes.binary[, .SD, .SDcols = setdiff(colnames(lst.checked_karyotypes.binary), c("MRN", "Normal", "accession_no", "clone.count", "chromosomal_count",  "ploidy", "full_deletion", "total_full_deletions", "MK_count", "MK", "CBF", "APL", "ck_3_or_more", "Incomplete.NA.NA", cbf.vec, apl.vec, who.vec))], 1, function(x) x <- sum(x, na.rm = T)))
  #  }
  #  lst.checked_karyotypes.binary$complex_BINARY  = ifelse(lst.checked_karyotypes.binary$abnormality_total > 2, 1, 0)
  #  lst.checked_karyotypes.binary$CBF <- unlist(apply(lst.checked_karyotypes.binary[, ..cbf.vec], 1, function(x) x <- ifelse(any(x), 1, 0)))
  #  lst.checked_karyotypes.binary$APL <- unlist(apply(lst.checked_karyotypes.binary[, ..apl.vec], 1, function(x) x <- ifelse(any(x), 1, 0)))
  #  lst.checked_karyotypes.binary$t9.11<- unlist(apply(lst.checked_karyotypes.binary[, ..t9.11.cols], 1, function(x) x <- ifelse(any(x), 1, 0)))
  #  lst.checked_karyotypes.binary$congenital <- unlist(apply(lst.checked_karyotypes.binary[, ..vec.congential.cols], 1, function(x) x <- ifelse(any(x), 1, 0)))
  
  
  # Accounting for NAs
  #  lst.checked_karyotypes.binary$MK = ifelse(is.na(lst.checked_karyotypes.binary$MK), 0, lst.checked_karyotypes.binary$MK)
  #  lst.checked_karyotypes.binary$complex_BINARY = ifelse(is.na(lst.checked_karyotypes.binary$complex_BINARY), 0, lst.checked_karyotypes.binary$complex_BINARY)
  #  lst.checked_karyotypes.binary$CBF = ifelse(is.na(lst.checked_karyotypes.binary$CBF), 0 , lst.checked_karyotypes.binary$CBF)
  #  lst.checked_karyotypes.binary$APL  = ifelse(is.na( lst.checked_karyotypes.binary$APL), 0,  lst.checked_karyotypes.binary$APL)
  #  lst.checked_karyotypes.binary$t9.11 = ifelse(is.na( lst.checked_karyotypes.binary$t9.11), 0,  lst.checked_karyotypes.binary$t9.11)
  #  lst.checked_karyotypes.binary$congenital = ifelse(is.na( lst.checked_karyotypes.binary$congenital), 0,  lst.checked_karyotypes.binary$congenital)
  
  # Labels
  #  lst.checked_karyotypes.binary[, cytogenetics := 'OND']
  #  if ("Normal" %in% names(lst.checked_karyotypes.binary)) {
  #    lst.checked_karyotypes.binary[Normal == T & abnormality_total == 0, cytogenetics := 'Normal']
  #  }
  #  lst.checked_karyotypes.binary[complex_BINARY == 1, cytogenetics := 'Complex']
  #  lst.checked_karyotypes.binary[t9.11 == 1, cytogenetics := 't9.11']
  #  lst.checked_karyotypes.binary[CBF == 1, cytogenetics := 'CBF']
  #  lst.checked_karyotypes.binary[APL == 1, cytogenetics := 'APL']
  #  lst.checked_karyotypes.binary[MK == 1, cytogenetics := 'Monosomal']
  
  #  lst.checked_karyotypes.binary$cytogenetics = factor( lst.checked_karyotypes.binary$cytogenetics , levels = c("CBF", "APL", "Monosomal", "Complex", "t9.11", "Normal", "OND"), labels = c("CBF", "APL", "Monosomal", "Complex", "t9.11", "Normal", "OND"))
  #  lst.checked_karyotypes.binary[!(cytogenetics %in% c("CBF", "APL")) & cell_total < 20, cytogenetics := factor("Unknown")]
  
  setnames(lst.checked_karyotypes.binary, "accession_no", "Accession.No")
  setcolorder(lst.checked_karyotypes.binary, c("MRN", "Accession.No", "chromosomal_count",  "ploidy", "clone.count", "cp"))
  
  return(lst.checked_karyotypes.binary)
  
}

#' @import data.table
NULL

#' Flags Incomplete Karyotypes
#'
#' Assigns feature type to each karyotype abnormality
#'
#' @param dt.path.cytoDmp.dx Table with parsed cytogenetics,  DMP panel, and cell count totals
#' @param dt.path.cyto.dx Table with report text
#' @param col.report.text Report Text
#' @param col.accession.no Accession No
#' @param col.MRN MRN column
#'
#' @return incomplete.mrns
#' @examples incomplete.mrns = flagIncompleteKaryotypes(dt.path.cytoDmp.dx, dt.path.cyto.dx)
#'
#'

flagIncompleteKaryotypes <- function(dt.path.cytoDmp.dx, dt.path.cyto.dx, lbl.risk.cyto.revised, col.report.text = "nan", col.accession.no = "Accession.No", col.MRN = "MRN", cell_lim = 20, col.cyto = "cytogenetics") {
  
  # Incomplete, cell count < 20
  #  dt.path.cyto.dx[, cell_total := sapply(get(col.report.text), cell_count_total)]
  
  #  merge.cols = c(col.MRN, col.accession.no)
  #  dt.incompletes = merge(dt.path.cytoDmp.dx, dt.path.cyto.dx[, .SD, .SDcols = c(merge.cols, "cell_total")], by = c(merge.cols))
  
  # Cell Count Total less than 20, excludes CBF and APL
  dt.incompletes = dt.path.cytoDmp.dx[cell_total < cell_lim,][!(get(col.cyto) %in% c("CBF", "APL")),] # CBF/APL will default to Good
  incomplete.mrns = unlist(dt.incompletes[, get(col.MRN)], use.names = F)
  
  # Incomplete, Normal/OND/t9.11 without DMP and no NCCN abnormalities
  abnormality.cols = cytoAbnormalities.vec.nonIntermediate(dt.path.cytoDmp.dx, lbl.risk.cyto.revised)
  
  dt.path.cytoDmp.dx[, nccn.abnormality_count := rowSums(dt.path.cytoDmp.dx[, .SD, .SDcols = (abnormality.cols)])]
  
  incomplete.mrns = c(incomplete.mrns, unlist(dt.path.cytoDmp.dx[nccn.abnormality_count == 0 | is.na(nccn.abnormality_count),][get(col.cyto) %in% c("Normal", "OND", "t9.11") & is.na(TP53), get(col.MRN)], use.names = F))
  
  #  incomplete.mrns = c(incomplete.mrns, unlist(dt.path.cytoDmp.dx[rowSums(dt.path.cytoDmp.dx[, ..abnormality.cols]) == 0][cytogenetics %in% c("Normal", "OND", "t9.11") & is.na(TP53), get(col.MRN)], use.names = F))
  
  #   incomplete.mrns = c(incomplete.mrns, unlist(dt.path.cytoDmp.dx[((cytogenetics %in% c("Normal", "OND", "t9.11")) & is.na(TP53)),  get(col.MRN)], use.names = F))
  
  # Normal/OND/t9.11 karyotype  w/o TP53 mutation and missing FLT3.ITD
  incomplete.mrns = c(incomplete.mrns, unlist(dt.path.cytoDmp.dx[nccn.abnormality_count == 0 | is.na(nccn.abnormality_count),][TP53 == F & is.na(FLT3.ITD) & get(col.cyto) %in% c("Normal", "OND", "t9.11"), get(col.MRN)], use.names= F))
  
  incomplete.mrns = unique(incomplete.mrns)
  return(incomplete.mrns)
}


#' Flag Karyotypes with donor cells
#'
#' Boolean indicator if the modal karyotype indicates donor cells
#'
#' @param str.ModalKaryotype Raw string with modal karyotype. Can be extracted using the `extractModalKaryotype` function
#'
#' @return TRUE/FALSE
#' @export

flagDonorKaryotype <- function(str.ModalKaryotype) {
  flag = grepl('\\/\\/', str.ModalKaryotype)
  return(flag)
}


#' Flags Incomplete Karyotypes (2022)
#'
#' Assigns feature type to each karyotype abnormality, ELN 2022
#'
#' @param dt.path.cytoDmp.dx Table with parsed cytogenetics,  DMP panel, and cell count totals
#' @param dt.path.cyto.dx Table with report text
#' @param col.report.text Report Text
#' @param col.accession.no Accession No
#' @param col.MRN MRN column
#'
#' @return incomplete.mrns
#' @examples incomplete.mrns = flagIncompleteKaryotypes.2022(dt.path.cytoDmp.dx, dt.path.cyto.dx)
#'

flagIncompleteKaryotypes.2022 <- function(dt.path.cytoDmp.dx, dt.path.cyto.dx, col.report.text = "nan", col.accession.no = "Accession.No", col.MRN = "MRN", cell_lim = 20, col.cyto = "cytogenetics") {
  
  # Cell Count Total less than 20, excludes CBF and APL
  dt.incompletes = dt.path.cytoDmp.dx[cell_total < cell_lim,][!(get(col.cyto) %in% c("CBF", "APL")),] # CBF/APL will default to Good
  incomplete.mrns = unlist(dt.incompletes[, get(col.MRN)], use.names = F)
  
  # Incomplete, Normal/OND/t9.11 without DMP and no NCCN abnormalities
  abnormality.cols = cytoAbnormalities.vec.nonIntermediate(dt.path.cytoDmp.dx)
  
  dt.path.cytoDmp.dx[, nccn.abnormality_count := rowSums(dt.path.cytoDmp.dx[, .SD, .SDcols = (abnormality.cols)])]
  
  incomplete.mrns = c(incomplete.mrns, unlist(dt.path.cytoDmp.dx[nccn.abnormality_count == 0 | is.na(nccn.abnormality_count),][get(col.cyto) %in% c("Normal", "OND", "t9.11") & is.na(TP53), get(col.MRN)], use.names = F))
  
  # Normal/OND/t9.11 karyotype  w/o TP53 mutation and missing FLT3.ITD
  incomplete.mrns = c(incomplete.mrns, unlist(dt.path.cytoDmp.dx[nccn.abnormality_count == 0 | is.na(nccn.abnormality_count),][TP53 == F  & get(col.cyto) %in% c("Normal", "OND", "t9.11")  & is.na(FLT3.ITD), get(col.MRN)], use.names= F))
  
  incomplete.mrns = unique(incomplete.mrns)
  return(incomplete.mrns)
}


#' Flags Mistyped Karyotypes
#'
#' Assigns feature type to each karyotype abnormality
#'
#' @param str.text Extracted Modal Karyotype from `extractModalKaryotype`
#'
#'
#' @return
#' @examples ISCN_flag = ISCN.check(str.text)
#' @export
#'

ISCN.check <- function(str.text) {
  num.check = as.integer(!grepl('\\[(cp)?(\\d+)](\\.|\\//|\\/)?$',  str.text))  + as.integer(grepl('[0-9]+[X|Y]',  str.text)  )
  if(num.check > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Flags Composite Karyotypes
#'
#' `TRUE` if karyotype is a composite
#'
#' @param str.text Extracted Modal Karyotype or cell count for the clone
#'
#'
#' @return
#' @examples cp.flag = composite.flag(str.text)
#' @export
#'

composite.flag <- function(str.text) {
  cp.flag = grepl('cp(\\d+)',  str.text)
  return(cp.flag)
}

#' @import data.table
NULL

#' Karyotype Indicator
#'
#' Logic, indicates whether a report string contains a karyotype
#'
#' @param text.example Pathology Report String
#'
#' @return TRUE or FALSE
#' @examples indicator = text.example(dt.path.orig[1, list(nan)])
#' @examples karyotype.isolated = text.example(dt.path.orig[1, list(nan)])
#' @export
#'


karyotype_included <- function(text.example) {
  
  indicator <- ifelse((grepl("Modal karyotype[:]", text.example) | grepl("Modal karyotype[(]s[)][:]", text.example)), T, F)
  text.example = text.example[grepl("Modal karyotype[:]|Modal karyotype[(]s[)][:]", text.example)]
  if (length(text.example) == 0) {
    return(F)
    break
  } else {
    if (grepl("Modal karyotype[:]", text.example)) {
      text.example = gsub(".*Modal karyotype[:]", "", text.example)  #Isolate Karyotype, left
    } else {
      text.example = gsub(".*Modal karyotype[(]s[)][:]", "", text.example)  #Isolate Karyotype, left
    }
    
    stop.words = c("Cytoscan.*|ANALYSIS.*|FISH.*|nuc.*|TEST.*|DISCLAIMER.*|DIAGNOSTIC.*|TestPerformed.*|METAPHASE.*|e-mailed.*|BONE.*|Note.*|Interim.*|PREPARATION.*|Probe.*|yes.*|ABNORMAL.*|not.*|\\*.*|random.*|Int(.)?.*")
    text.example = gsub(stop.words, "", text.example)
    
    #  text.example = gsub("Cytoscan.*|\\*.*|ANALYSIS.*|FISH.*|nuc.*|TEST.*|DISCLAIMER.*|DIAGNOSTIC.*|e-mailed.*|BONE.*|Note.*|Interim.*|PREPARATION.*|Probe.*|yes.*|ABNORMAL.*|not.*|Int(.)?.*", "", text.example)
    
    text.example = gsub('^/', '', text.example) #remove leading blank
    text.example = gsub('\\/', '', text.example) #remove slash
    text.example = gsub("\\n", "", text.example, perl = T) # Remove Carriage Returns
    text.example = gsub("\\\\n", "", text.example)
    text.example = gsub("idemdel", "idem,del", text.example) # Remove Typos
    text.example = gsub("\\s", "", text.example) # Remove White Spaces
    text.example = gsub("\\?", "", text.example) # Remove Question Mark
    
    
    text.example = unlist(strsplit(text.example, "]")) # Split on square bracket , unlist
    text.example = paste0(text.example, "]") # Add back bracket
    text.example = setdiff(text.example, "]")
    text.example = setdiff(text.example, ".]")
    if (length(text.example) == 0) { # For bad reports :(
      return(F)
    } else {
      return(T)
    }
  }
}

#' Extract Karyotype
#'
#' Extracts Raw Karyotype from Report Text (No parsing/cell counts)
#'
#' @param text.example Pathology Report String
#'
#' @return karyotype.raw
#' @examples karyotype.raw = extractModalKaryotype("Pathology report--Modal karyotype: 45,XX,-5[20]")
#' @examples karyotype.raw = "45,XX,-5[20]"
#' @export
#'


extractModalKaryotype <- function(text.example) {
  text.example = as.character(text.example)
  text.example = cleanPath.str(path.str = text.example)
  text.example = text.example[grepl("Modal karyotype[:]|Modal karyotype[(]s[)][:]",
                                    text.example)]
  if (length(text.example) == 0) {
    return(NA)
    break
  }
  if (grepl("Modal karyotype[:]", text.example)) {
    text.example = gsub(".*Modal karyotype[:]", "", text.example)
  }
  else {
    text.example = gsub(".*Modal karyotype[(]s[)][:]", "",
                        text.example)
  }
  #  text.example = gsub("Cytoscan.*|ANALYSIS.*|FISH.*|nuc.*|TEST.*|DISCLAIMER.*|DIAGNOSTIC.*|TestPerformed.*|METAPHASE.*|e-mailed.*|BONE.*|Note.*|Interim.*|PREPARATION.*|Probe.*|yes.*|ABNORMAL.*|not.*|\\*.*|random.*|Int(.)?.*",
  #                      "", text.example, ignore.case = T)
  stop.words = c("Cytoscan.*|ANALYSIS.*|FISH.*|nuc.*|TEST.*|DISCLAIMER.*|DIAGNOSTIC.*|TestPerformed.*|METAPHASE.*|e-mailed.*|BONE.*|Note.*|Interim.*|PREPARATION.*|Probe.*|yes.*|ABNORMAL.*|not.*|\\*.*|random.*|Int(.)?.*|AMENDED.*")
  text.example = gsub(stop.words, "", text.example)
  
  text.example = gsub("^/", "", text.example)
  #  text.example = gsub("\\/", "", text.example)
  text.example = gsub("\\n", "", text.example, perl = T)
  text.example = gsub("\\\\n", "", text.example)
  text.example = gsub("idemdel", "idem,del", text.example)
  text.example = gsub("\\s", "", text.example)
  text.example = gsub("\\?", "", text.example)
  #  text.example = unlist(strsplit(text.example, "]"))
  #  text.example = paste0(text.example, "]")
  #  text.example = setdiff(text.example, "]")
  #  text.example = setdiff(text.example, ".]")
  if (length(text.example) == 0) {
    return(NA)
  }
  else {
    return(text.example)
  }
}


#' Isolate Karyotype
#'
#' Extracts Karyotype from the Report Text
#'
#' @param text.example Pathology Report String
#'
#' @return karyotype.isolated
#' @examples karyotype.isolated = text.example(dt.path.orig[1, list(nan)])
#' @export
#'


karyotype_parser <- function(text.example)  {
  text.example = as.character(text.example)
  text.example = cleanPath.str (path.str = text.example)
  text.example = text.example[grepl("Modal karyotype[:]|Modal karyotype[(]s[)][:]", text.example)]
  if (length(text.example) == 0) {
    return(NA)
    break
  }
  
  if (grepl("Modal karyotype[:]", text.example)) {
    text.example = gsub(".*Modal karyotype[:]", "", text.example)  #Isolate Karyotype, left
  } else {
    text.example = gsub(".*Modal karyotype[(]s[)][:]", "", text.example)  #Isolate Karyotype, left
  }
  
  #  text.example = gsub("Cytoscan.*|ANALYSIS.*|FISH.*|nuc.*|TEST.*|DISCLAIMER.*|DIAGNOSTIC.*|TestPerformed.*|METAPHASE.*|e-mailed.*|BONE.*|Note.*|Interim.*|PREPARATION.*|Probe.*|yes.*|ABNORMAL.*|not.*|random.*|\\*.*|Int(.)?.*", "", text.example, ignore.case = T)
  stop.words = c("Cytoscan.*|ANALYSIS.*|FISH.*|nuc.*|TEST.*|DISCLAIMER.*|DIAGNOSTIC.*|TestPerformed.*|METAPHASE.*|e-mailed.*|BONE.*|Note.*|Interim.*|PREPARATION.*|Probe.*|yes.*|ABNORMAL.*|not.*|\\*.*|random.*|AMENDED.*|Int(.)?.*")
  text.example = gsub(stop.words, "", text.example,  ignore.case = T)
  
  text.example = gsub('^/', '', text.example) #remove leading blank
  text.example = gsub('\\/', '', text.example) #remove slash
  text.example = gsub("\\n", "", text.example, perl = T) # Remove Carriage Returns
  text.example = gsub("\\\\n", "", text.example)
  text.example = gsub("idemdel", "idem,del", text.example) # Remove Typos
  text.example = gsub("\\s", "", text.example) # Remove White Spaces
  text.example = gsub("\\?", "", text.example) # Remove Question Mark
  
  
  text.example = unlist(strsplit(text.example, "]")) # Split on square bracket , unlist
  text.example = paste0(text.example, "]") # Add back bracket
  text.example = setdiff(text.example, "]")
  text.example = setdiff(text.example, ".]")
  if (length(text.example) == 0) { # For bad reports :(
    return(NA)
  } else {
    
    lst.text = split(text.example, 1:length(text.example)) # split karyotype into list
    
    cell.counts = gsub("\\[|\\]", "", regmatches(lst.text, gregexpr("\\[.*?\\]", lst.text))) # Extract Cell Counts
    # cell.counts = gsub("cp", "", cell.counts) # Swap out CP
    names(cell.counts) = rep("Cell Counts", length(cell.counts)) # Name Cell Counts Vector
    
    lst.clones = split(cell.counts, 1:length(text.example)) # Create list w cell counts
    names(lst.clones) = split(paste0("Clone ", 1:length(text.example)), 1:length(text.example)) # Name Clone List
    
    text.example = gsub("\\[[^()]*\\]", "", text.example) # Remove Cell Count
    
    dt.example  = data.table(karyotype = text.example) # put on separate rows
    dt.example[, clone := paste("Clone", 1:nrow(dt.example))]
    dt.example[, idem := grepl("idem|idem\n|idem?x?[0-9]", karyotype)] #idem indicator
    dt.example[, sdl := grepl("sdl|sdl\n|sdl?[0-9]|sl|sl\n", karyotype)] #sideline indicator
    dt.example[, idem_expanded := ifelse(idem == T, gsub(".*XX,|.*XY,", " ", karyotype)[which(idem == F)][1], NA)]
    dt.example[, parent := ifelse(idem == TRUE, dt.example$clone[dt.example$idem == F][1], NA)] #Parent for idem
    dt.example[, parent := ifelse(sdl == T, shift(clone, n = 1), parent)] # Parent for sdl
    dt.example[, ploidy := sapply(karyotype, function(x)  x = regmatches(x, regexpr("<[^>]+>", x))[1])] # Add ploidy info
    dt.example[is.na(ploidy), ploidy:='<2n>'] # Defaulted to 2N
    
    #  dt.example[, cell.counts := cell.counts] # Add cell counts
    #  dt.example[, cp := grepl("cp", cell.counts)] # Add composite info
    
    text.example = strsplit(text.example, ",")
    
    
    # Add Cell Counts in
    for (i in 1:length(text.example)) {
      for (j in length(text.example[[i]])) {
        text.example[[i]][j+1] = paste0("Cell Counts = ", cell.counts[i])
        text.example[[i]][j+2] = paste0("Parent = ", dt.example$parent[i])
        text.example[[i]][j+3] = paste0("Ploidy = ", dt.example$ploidy[i])
        text.example[[i]][j+4] = paste0("Composite = ", grepl("cp", cell.counts[i]))
      }
    }
    
    
    names(text.example) = paste0("Clone ", 1:length(text.example))
    karyotype.isolated = text.example
    return(karyotype.isolated) }
}


#' Isolate Karyotype (Not in a Pathology Report)
#'
#' Extracts Karyotype from a character string
#'
#' @param text.example Pathology Report String
#'
#' @return karyotype.isolated
#' @examples karyotype.isolated = karyotype_parser.non_report("46,XX[18]/47,XX,+4[2]")
#' @export
#'


karyotype_parser.non_report <- function(text.example)  {
  text.example = as.character(text.example)
  text.example = paste0("Modal karyotype:", text.example)
  karyotype.isolated = karyotype_parser(text.example)
  return(karyotype.isolated)
}


#' Flags Karyotype Features
#'
#' Assigns feature type to each karyotype abnormality
#'
#' @param karyotype.isolated Isolated karyotype from `karyotype_parser`
#'
#' @return lst.abnormalities
#' @examples `lst.abnormalities = flag_function(karyotype.isolated)`
#' @export
#'


flag_function <- function(text.example) {
  
  lst.abnormalities = text.example
  for (i in 1:length(text.example)) { #iterating for each clone
    for (j in 1:length(text.example[[i]])) { #Iterating for each feature
      #     feature = text.example[[i]][j]
      feature = text.example[[i]][j]
      if(grepl("Ploidy", feature))  {
        lst.feature = as.list(c("ploidy" =  unlist(strsplit(feature, " = "))[2]))
        lst.abnormalities[[i]][j] = list(lst.feature)
        next
      } else {
        feature = gsub("<[^>]+>", "", feature, perl = T)  # Remove Text in angle brackets
      }
      if (grepl("Cell Counts", feature)) { # cell counts
        lst.feature = as.list(c("cell_count" = unlist(strsplit(feature, " = "))[2]))
      } else  if (grepl("Composite", feature)) {
        lst.feature = as.list(c("cp" = unlist(strsplit(feature, " = "))[2]))
      } else if (grepl("Parent", feature)) { #Identify Parent
        parent = unlist(strsplit(feature, " = "))[2]
        if (grepl("NA", parent)) {
          lst.feature = as.list(c("parent_clone" = "No Parent"))
        } else {
          lst.feature = as.list(c("parent_clone" = parent))
        }
      } else if (grepl("idem|idem\n|idem?x?[0-9]|sdl|sdl\n|sdl?[0-9]|sl|sl\n", feature)) { #idem, sdl, sdl
        lst.feature = as.list(c("subclone" = feature))
      } else if (!(grepl("[a-z]+|^-|\\+|[A-Z]+", feature))) { # Not an abnormality or sex chromosomes (Chrom Count)
        lst.feature = as.list(c("chromosomal_count" = feature))
      } else if (grepl("c", feature) & !(grepl("inc", feature))) { # Congential abnormalities
        chrom.cong = regmatches(feature, regexpr("[0-9]+", feature))[1]
        chrom.region = gsub("\\(|\\)", "", regmatches(feature, gregexpr("\\(+[p|q]", feature)))
        if (chrom.region == "character0") {
          chrom.region <- "FULL"
        }
        lst.feature = as.list(c("abnormality" = "congenital", "chrom1" = chrom.cong, "region1" = chrom.region))
      } else if (grepl("dmin", feature)) { # Minute
        lst.feature = as.list(c("abnormality" = "double-minute", "count" = gsub("dmin", "\\1", feature)))
      } else if (grepl("\\d+\\-\\d+", feature) & !(grepl("mar", feature))) { # Chromosomal range count
        lst.feature = as.list(c("chromosomal_count" = feature))
      } else if (grepl("XY|XX|Y|X", feature) & !(grepl("[+|-][X|Y]", feature)) & !(grepl("t\\(", feature))) { # Sex Chromosomes, excluding deletions/additions, no sex chromosome translocations
        lst.feature = as.list(c("sex_chrom" = feature))
      } else if (grepl("del", feature) & !grepl("der", feature)) { # Partial Deletion
        del.chrom = gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", gsub("([^()]+\\(([0-9]+|X|Y)\\))\\([^)]+\\)", "\\1", feature), perl = T)
        del.chrom =  unlist(regmatches(del.chrom, gregexpr("([0-9]+|X|Y)", del.chrom)))
        #  del.region = gsub("\\(|\\)", "", regmatches(feature, gregexpr("\\(+[p|q]", feature)))
        del.region = unlist(regmatches(feature, gregexpr("[p|q]", feature)))
        del.region = del.region[1] # Interstitial Deletions
        lst.feature = as.list(c("abnormality" = "deletion", "chrom1" = del.chrom, "region1" = del.region))
      } else if (grepl("inv", feature)) { # Inversion
        inv.chrom = gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", gsub("([^()]+\\([0-9]+\\))\\([^)]+\\)", "\\1", feature), perl = T)
        inv.region = unlist(regmatches(feature, gregexpr("[p|q]", feature)))
        lst.feature = as.list(c("abnormality" = "inversion", "chrom1" = inv.chrom, "region" = inv.region))
      } else if (grepl("^[-]", feature) & !grepl("add|del|mar|r", feature)) { # Full Deletion
        lst.feature = as.list(c("abnormality" = "deletion", "chrom1"= gsub("[-]", "",feature), "region1" = "FULL"))
      } else if (grepl("add", feature) & !grepl("der", feature)) { # Partial Addition
        add.chrom = gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", gsub("([^()]+\\([0-9]+\\))\\([^)]+\\)", "\\1", feature), perl = T)
        add.region = gsub("\\(|\\)", "", regmatches(feature, gregexpr("\\(+[p|q]", feature)))
        lst.feature = as.list(c("abnormality" = "addition", "chrom1" = add.chrom, "region1" = add.region))
      } else if (grepl("mar?x?[0-9]|?[0-9]mar|mar|?\\+mar", feature)) { # Marker Chromosomes
        #        lst.feature = as.list(c("abnormality" = "marker", "count"= gsub("mar", "\\1", feature)))
        marker.count = regmatches(feature, gregexpr("[0-9]+", feature))
        marker.count = ifelse(length(unlist(marker.count)) == 0, 1, marker.count)
        lst.feature = as.list(c("abnormality" = "marker", "count"= marker.count))
      } else if (grepl("der", feature)) { # Derivative Chromosome
        der.chrom = gsub("\\(|\\)", "", regmatches(feature, gregexpr("der\\([0-9]+\\)", feature)))
        der.chrom = gsub(".*der", "", der.chrom)
        if (grepl("t\\(", feature)) { # If Derivative has Translocation
          translocation.part1 = gsub("", "", regmatches(feature, gregexpr("t\\([0-9|X]+\\;[0-9|X|Y]+\\)", feature))) # Extract translocation
          translocation.part1 = unlist(regmatches(translocation.part1, gregexpr("\\d+|X", translocation.part1)))
          translocation.part2 = gsub(".*t\\([0-9|X|Y]+\\;[0-9|X|Y]+\\)",  "\\1", feature, perl = T) # Extract p/q arms
          translocation.part2 = unlist(regmatches(translocation.part2, gregexpr("[p|q][0-9]+", translocation.part2)))
          lst.feature = as.list(c("abnormality" = "derivative", "chrom1" = der.chrom, "region1" = translocation.part1, "region2" = translocation.part2))
        } else if (der.chrom == "character0")  {
          der.part1 = gsub("", "", regmatches(feature, gregexpr("der\\([0-9]+\\;[0-9]+\\)", feature))) # Extract translocation
          der.part1 = unlist(regmatches(der.part1, gregexpr("\\d+", der.part1)))
          der.part2 = gsub(".*der\\([0-9]+\\;[0-9]+\\)",  "\\1", feature, perl = T) # Extract p/q arms
          der.part2 = unlist(regmatches(der.part2, gregexpr("[p|q][0-9]+", der.part2)))
          lst.feature = as.list(c("abnormality" = "derivative", "chrom1" = "UNK", "region1" = der.part1, "region2" = der.part2))
        } else {
          lst.feature = as.list(c("abnormality" = "derivative", "chrom1"= der.chrom))
        }
        if (grepl("hsr", feature)) { #HSR
          lst.feature[["hsr"]] = "hsr"
        }
      } else if (grepl("t\\([0-9|X|Y]+\\;[0-9|X|Y]+\\;[0-9|X|Y]+\\)", feature)) {  #Translocation, 3 chromosomes
        translocation.part1 = gsub("", "", regmatches(feature, gregexpr("t\\([0-9|X|Y]+(\\;|\\:)[0-9|X|Y]+\\;[0-9|X|Y]+\\)", feature))) # Extract translocation
        translocation.part1 =unlist(regmatches(translocation.part1, gregexpr("\\d+|X|Y", translocation.part1)))
        translocation.part2 = gsub(".*t\\([0-9|X|Y]+\\;[0-9|X|Y]+\\;[0-9|X|Y]+\\)",  "\\1", feature, perl = T) # Extract p/q arms
        translocation.part2 = unlist(regmatches(translocation.part2, gregexpr("[p|q][0-9]+", translocation.part2)))
        lst.feature = as.list(c("abnormality" = "translocation", "chrom" = translocation.part1, "region" = translocation.part2))
      } else if (grepl("t\\([0-9|X|Y]+\\;[0-9|X|Y]+\\;[0-9|X|Y]+\\;[0-9|X|Y]+\\)", feature)) {  #Translocation, 4 chromosomes :(
        translocation.part1 = gsub("", "", regmatches(feature, gregexpr("t\\([0-9|X|Y]+\\;[0-9|X|Y]+\\;[0-9|X|Y]+\\;[0-9|X|Y]+\\)", feature))) # Extract translocation
        translocation.part1 =unlist(regmatches(translocation.part1, gregexpr("\\d+|X|Y", translocation.part1)))
        translocation.part2 = gsub(".*t\\([0-9|X|Y]+\\;[0-9|X|Y]+\\;[0-9|X|Y]+\\;[0-9|X|Y]+\\)",  "\\1", feature, perl = T) # Extract p/q arms
        translocation.part2 = unlist(regmatches(translocation.part2, gregexpr("[p|q][0-9]+", translocation.part2)))
        lst.feature = as.list(c("abnormality" = "translocation", "chrom" = translocation.part1, "region" = translocation.part2))
      } else if (grepl("t\\(", feature)) { # Translocation, 2 chromosomes
        translocation.part1 = gsub("", "", regmatches(feature, gregexpr("t\\([0-9|X|Y]+(\\;|\\:)[0-9|X]+\\)", feature, perl = T))) # Extract translocation (including X chromosomes)
        translocation.part1 = unlist(regmatches(translocation.part1, gregexpr("\\d+|X|Y", translocation.part1)))
        translocation.part2 = gsub(".*t\\([0-9|X|Y]+\\;[0-9|X|Y]+\\)",  "\\1", feature, perl = T) # Extract p/q arms
        translocation.part2 = unlist(regmatches(translocation.part2, gregexpr("[p|q][0-9]+", translocation.part2)))
        lst.feature = as.list(c("abnormality" = "translocation", "chrom" = translocation.part1, "region" = translocation.part2))
      } else if (grepl("^[+][0-9]+", feature)) { # Full Addition
        lst.feature = as.list(c("abnormality" = "addition", "chrom1"= gsub("[+]", "",feature), "region1" = "FULL"))
      } else if (grepl("(?<!i)dic", feature, perl = T)) { # dicentric, not isodicentric
        dic.part1 = gsub("", "", regmatches(feature, gregexpr("dic\\([0-9]+\\;[0-9]+\\)", feature))) # Extract translocation
        dic.part1 = unlist(regmatches(dic.part1, gregexpr("\\d+", dic.part1)))
        dic.part2 = gsub(".*dic\\([0-9]+\\;[0-9]+\\)",  "\\1", feature, perl = T) # Extract p/q arms
        dic.part2 = unlist(regmatches(dic.part2, gregexpr("[p|q][0-9]+", dic.part2)))
        lst.feature = as.list(c("abnormality" = "dicentric", "chrom" = dic.part1, "region" = dic.part2))
      }  else if (grepl("idic", feature)) { # isodicentric
        dic.part1 = gsub("", "", regmatches(feature, gregexpr("idic\\([0-9]+\\)", feature))) # Extract first part
        dic.part1 = unlist(regmatches(dic.part1, gregexpr("\\d+", dic.part1)))
        dic.part2 = gsub(".*idic\\([0-9]+\\)",  "\\1", feature, perl = T) # Extract p/q arms
        dic.part2 = unlist(regmatches(dic.part2, gregexpr("[p|q][0-9]+", dic.part2)))
        lst.feature = as.list(c("abnormality" = "isodicentric", "chrom" = dic.part1, "region" = dic.part2))
      } else if (grepl("Parent|Child", feature)) {
        lst.feature = as.list(c("Clone Type" = feature))
      } else if (grepl("inc", feature)) { # Incomplete
        lst.feature = as.list(c("abnormality" = "Incomplete"))
      } else if (grepl("dup", feature)){ # Duplication
        dup.part1 = gsub("", "", regmatches(feature, gregexpr("dup\\([0-9]+\\)", feature)))
        dup.part1 = unlist(regmatches(dup.part1, gregexpr("\\d+", dup.part1)))
        dup.part2 = gsub(".*dup\\([0-9]+\\)",  "\\1", feature, perl = T) # Extract p/q arms
        dup.part2 = unlist(regmatches(dup.part2, gregexpr("[p|q][0-9]+", dup.part2)))
        lst.feature = as.list(c("abnormality" = "duplication", "chrom1" = dup.part1, "region" = dup.part2))
      } else if (grepl("i\\([0-9]+?[p|q]\\)", feature)) { # Isochromosome
        i.chrom = gsub("", "", regmatches(feature, gregexpr("i\\([0-9]+?[p|q]\\)", feature)))
        i.chrom = unlist(regmatches(i.chrom, gregexpr("\\d+", i.chrom)))
        i.region = gsub(".*i\\([0-9]+?[p|q]\\)",  "\\1", feature, perl = T) # Extract p/q arms
        i.region = unlist(regmatches(i.region, gregexpr("[p|q][0-9]+", i.region)))
        lst.feature = as.list(c("abnormality" = "isocentric", chrom1 = i.chrom, region1 = i.region))
      } else if (grepl("^[+][X|Y]+", feature)) { # Full Addition (Sex Chromosomes)
        lst.feature = as.list(c("abnormality" = "addition", "chrom1"= gsub("[+]", "",feature), "region1" = "FULL"))
      } else if (grepl("^[-][X|Y]+", feature)) { # Full Deletion (Sex Chromosomes)
        lst.feature = as.list(c("abnormality" = "deletion", "chrom1"= gsub("[-]", "",feature), "region1" = "FULL"))
      } else if (grepl("ins\\([0-9]+\\;?", feature)) { # Insertion w/unknown
        ins.chrom = gsub("", "", regmatches(feature, gregexpr("ins\\([0-9]+\\;", feature)))
        ins.chrom = unlist(regmatches(ins.chrom, gregexpr("\\d+", ins.chrom)))
        ins.region = gsub(".*ins\\([0-9]+?;\\)",  "\\1", feature, perl = T) # Extract p/q arms
        ins.region = unlist(regmatches(ins.region, gregexpr("[p|q][0-9]+", ins.region)))
        lst.feature = as.list(c("abnormality" = "insertion", chrom1 = ins.chrom, region1 = ins.region))
      } else if (grepl("?(+|-|~)r", feature)) { # Ring Chromosome
        ring.chrom = gsub("", "", regmatches(feature, gregexpr("r\\([0-9]+\\)", feature)))
        ring.chrom = unlist(regmatches(ring.chrom, gregexpr("\\d+", ring.chrom)))
        ring.region = gsub("\\(|\\)", "", regmatches(feature, gregexpr("\\(+[p|q]", feature)))
        lst.feature = as.list(c("abnormality" = "ring", chrom1 = ring.chrom))
      } else  {
        lst.feature = as.list(c("abnormality" = "UNKNOWN", feature = feature))
      }
      if (is.null(lst.feature)) {
        lst.abnormalities[[i]][j] = NULL
      } else {
        lst.abnormalities[[i]][j] = list(lst.feature)
      }
    }
  }
  return(lst.abnormalities)
}


#' Tabulates Karyotype Features
#'
#' Converts flagged list to a feature table. Allows for filtering of clones with low cell counts.
#'
#' @param lst.abnormalities List of flagged abnormalities from `flag_function`
#' @param MRN Individual MRN
#' @param accession.no Individual Accession NO
#' @param cell_count_lim The cell count limit for aggregation. Default is 1
#'
#' @return dt.clones
#' @examples dt.clones = rbindlist.abnormalities(lst.abnormalities, MRN = "#000000", accession.no = "BBBCS-1216")
#' @export
#'

rbindlist.abnormalities <- function(lst.abnormalities, MRN, accession.no, cell_count_lim = 1) {
  
  
  lst.dt.abnormalities= list()
  for (i in 1:length(lst.abnormalities)) {
    dt.abnormalities = data.table(1)[,`:=`(c("feature", "value", "feature_no"),NA)][,V1:=NULL][.0]
    for (j in 1:length(lst.abnormalities[[i]])) {
      # lst.abnormalities[[i]]
      vec.abnormal = unlist(lst.abnormalities[[i]][j], use.names = T)
      df.abnormalities = data.frame(feature = names(vec.abnormal), value = vec.abnormal)
      dt.abnormalities.sub = data.table(df.abnormalities)
      dt.abnormalities.sub[, feature_no := j]
      dt.abnormalities = rbind(dt.abnormalities.sub, dt.abnormalities)
    }
    
    dt.features = dt.abnormalities[!(feature %in% c("chromosomal_count", "sex_chrom", "parent_clone", "cell_count", 'ploidy', "clone_no",  "cp")),]
    if (nrow(dt.features) == 0) {
      dt.features = data.table()
      dt.features[, chromosomal_count := ifelse(length(dt.abnormalities[feature == "chromosomal_count", value]) > 1, dt.abnormalities[feature == "chromosomal_count", value][-1], dt.abnormalities[feature == "chromosomal_count", value])]
      dt.features[, sex_chrom := ifelse(length(dt.abnormalities[feature == "sex_chrom", value]) > 1, dt.abnormalities[feature == "sex_chrom", value][-1], dt.abnormalities[feature == "sex_chrom", value])]
      dt.features[, parent_clone := dt.abnormalities[feature == "parent_clone", value]]
      dt.features[, clone_no := rep(names(lst.abnormalities[i]), nrow(dt.features))]
      dt.features[, cell_count := dt.abnormalities[feature == "cell_count", value]]
      dt.features[, ploidy := dt.abnormalities[feature == "ploidy", value]]
      dt.features[, cp := dt.abnormalities[feature == "cp", value]]
      lst.dt.abnormalities[[i]] = dt.features
      
    } else {
      dt.features = dcast(dt.features[, list(feature_no, feature, value)], feature_no ~ feature)
      
      dt.features[, chromosomal_count := ifelse(length(dt.abnormalities[feature == "chromosomal_count", value]) > 1, dt.abnormalities[feature == "chromosomal_count", value][-1], dt.abnormalities[feature == "chromosomal_count", value])]
      dt.features[, sex_chrom := ifelse(length(dt.abnormalities[feature == "sex_chrom", value]) > 1, dt.abnormalities[feature == "sex_chrom", value][-1], dt.abnormalities[feature == "sex_chrom", value])]
      
      # dt.features[, sex_chrom := dt.abnormalities[feature == "sex_chrom", value][-1]]
      dt.features[, parent_clone := dt.abnormalities[feature == "parent_clone", value]]
      dt.features[, clone_no := rep(names(lst.abnormalities[i]), nrow(dt.features))]
      dt.features[, cell_count := dt.abnormalities[feature == "cell_count", value]]
      dt.features[, ploidy := dt.abnormalities[feature == "ploidy", value]]
      dt.features[, cp := dt.abnormalities[feature == "cp", value]]
      
    }
    # if (nrow(dt.features[cell_count > 1, ]) == 0) { # Remove clones with one cell count
    #    next
    # } else {
    lst.dt.abnormalities[[i]] = dt.features
    #  }
  }
  
  dt.clones = rbindlist(lst.dt.abnormalities, fill = T)
  dt.clones[, MRN :=  MRN]
  dt.clones[, accession_no :=  accession.no]
  dt.clones.names = c("MRN", "accession_no", "clone_no", "chromosomal_count", "sex_chrom", "parent_clone", "cell_count", "ploidy", "cp", "abnormality", "chrom1" ,"region", "region11", "region12", "region21", "region22")
  dt.clones.names = intersect(dt.clones.names, names(dt.clones))
  setcolorder(dt.clones, dt.clones.names)
  #   dt.clones =  dt.clones[cell_count > 1, ]
  # accounting for cell count lim but is a parent
  dt.clones[, drop_count := ifelse(!(dt.clones$clone_no %in% dt.clones$parent_clone) & dt.clones$cell_count == cell_count_lim, T, F)]
  dt.clones = dt.clones[drop_count == F,]
  dt.clones[, drop_count:= NULL]
  # Remove NA columns from cell count 1
  dt.clones = dt.clones[, .SD, .SDcols = names(which(dt.clones[, colSums(is.na(dt.clones)) != nrow(dt.clones)]))]
  if (nrow(dt.clones) == 0) {
    vec.message = paste0("Cell Count = ", cell_count_lim)
    return(vec.message)
  } else {
    return(dt.clones)
  }
}


#' Aggregates Karyotype Abnormalities
#'
#' Converts feature table into summary table
#'
#' @param kary.table Features table from `rbindlist.abnormalities`
#'
#' @return kary.table.dcast
#' @examples tbl.dcast = aggregate_tbl(features.tbl)
#' @export
#'


aggregate_tbl <- function(kary.table) {
  if (!is.data.table(kary.table)){
    kary.table.dcast = kary.table
    return(kary.table.dcast)
  }
  if (!is.null(kary.table$subclone)) {
    kary.table.sub = kary.table[is.na(subclone) & !is.na(abnormality),]
  } else if (!is.null(kary.table$abnormality)) {
    if("region1" %in% names(kary.table)) {
      full.deletion.check.clones = kary.table[as.numeric(cell_count) < 3 & abnormality == 'deletion' & region1 == 'FULL', clone_no]
      if(length(full.deletion.check.clones) > 0) {
        full.deletion = kary.table[clone_no %in% full.deletion.check.clones, ]
        full.deletion = full.deletion[,.N, by = 'clone_no']
        full.deletion.insufficient = full.deletion[N == 1, clone_no]
        vec.parent.clone.subset = kary.table[clone_no %in% parent_clone] # accounting if the clone is a parent (clone should not be a parent to get knocked off)
        if(length(full.deletion.insufficient) > 0 & nrow(vec.parent.clone.subset) == 0){
          kary.table.sub = kary.table[!(clone_no %in% full.deletion.insufficient)]
          if(nrow(kary.table.sub) == 0) {
            vec.message = paste0("Cell Count = ", 2)
            return(vec.message)
          }
        } else {
          kary.table.sub = kary.table[!is.na(abnormality),]
        }
      } else {
        kary.table.sub = kary.table[!is.na(abnormality),]
      }
    }
    else {
      kary.table.sub = kary.table[!is.na(abnormality),]
    }
  } else {
    kary.table.sub = kary.table
  }
  if (is.null(kary.table.sub$abnormality)) {
    kary.table.sub[, col_name := "Normal"]
    kary.table.dcast = dcast(kary.table.sub[, list(accession_no, MRN, col_name)], MRN + accession_no ~ col_name, fun=length)
    if("chromosomal_count" %in% names(kary.table)) {
      kary.table.dcast[, chromosomal_count := toString(unique(kary.table[, (chromosomal_count)]))]
    }
  } else {
    if (is.null(kary.table.sub$count) & !(is.null(kary.table.sub$count1) & is.null(kary.table.sub$count2))) {
      kary.table.sub[, count := paste0(count1, "_", count2)]
    }
    if (c("chrom1") %in% colnames(kary.table.sub)) {
      kary.table.sub[, chrom_var := chrom1]
      kary.table.sub[(abnormality == "marker" | abnormality == "double-minute"), chrom_var := chrom1]
      kary.table.sub[(abnormality == "translocation" & !("chrom3" %in% names(kary.table.sub))), chrom_var := paste0(chrom1, ".", chrom2)]
      kary.table.sub[(abnormality == "translocation" & ("chrom3" %in% names(kary.table.sub))), chrom_var := paste0(chrom1, ".", chrom2, ".", chrom3)]
      kary.table.sub[(abnormality == "Incomplete" | abnormality == "UNKNOWN"), chrom_var := NA]
    } else {
      kary.table.sub[, chrom_var := NA]
      kary.table.sub[abnormality == "marker" | abnormality == "double-minute", chrom_var := count]
      kary.table.sub[abnormality == "translocation" & !("chrom3" %in% names(kary.table.sub)), chrom_var := paste0(chrom1, ".", chrom2)]
      kary.table.sub[abnormality == "translocation" & ("chrom3" %in% names(kary.table.sub)), chrom_var := paste0(chrom1, ".", chrom2, ".", chrom3)]
      kary.table.sub[(abnormality == "Incomplete" | abnormality == "UNKNOWN"), chrom_var := NA]
      
    }
    kary.table.sub[, region_var := ""]
    kary.table.sub[(abnormality == "translocation" | abnormality == "inversion" | abnormality == "duplication")  & length(intersect(names(kary.table.sub), c("region1", "region2"))) == 2, region_var := paste0(region1, ".", region2)]
    if (length(intersect("region3",names(kary.table.sub))) == 1)  {
      # kary.table.sub[(abnormality == "translocation")  & length(intersect(names(kary.table.sub), c("region1", "region2", "region3"))) == 3, region_var := paste0(region1, ".", region2, ".", region3)]
      kary.table.sub[(abnormality == "translocation")  & !is.na(region3), region_var := paste0(region1, ".", region2, ".", region3)]
    }
    kary.table.sub[abnormality == "derivative" & length(intersect(names(kary.table.sub), c("region1", "region2"))) == 2, region_var := paste0(region1, ".", region2)]
    kary.table.sub[abnormality == "derivative" & length(intersect(names(kary.table.sub), c("region11", "region12"))) == 2, region_var := paste0(region11, ".", region12)]
    kary.table.sub[abnormality == "Incomplete" | abnormality == "UNKNOWN", region_var := NA]
    kary.table.sub[abnormality == "marker", region_var := ""]
    kary.table.sub[(abnormality %in% c("deletion", "addition", "ring") & length(intersect(names(kary.table.sub), c("region1"))) == 1), region_var := region1]
    kary.table.sub[(length(intersect(names(kary.table.sub), c("region1"))) == 1) & region_var == "", region_var := region1]
    
    kary.table.sub[, col_name := paste(abnormality,chrom_var, region_var, sep = ".")]
    
    # Full Deletion check
    full.deletion.check.clones = kary.table.sub[as.numeric(cell_count) < 3 & grepl('deletion\\.\\d+\\.FULL', col_name), clone_no]
    if(length(full.deletion.check.clones) > 0) {
      full.deletion = kary.table.sub[clone_no %in% full.deletion.check.clones, ]
      full.deletion = full.deletion[,.N, by = 'clone_no']
      full.deletion.insufficient = full.deletion[N == 1, clone_no]
      vec.parent.clone.subset = kary.table.sub[clone_no %in% parent_clone] # accounting if the clone is a parent (clone should not be a parent to get knocked off)
      if(length(full.deletion.insufficient) > 0 & nrow(vec.parent.clone.subset) == 0){
        kary.table.sub = kary.table.sub[!(clone_no %in% full.deletion.insufficient)]
        if(nrow(kary.table.sub) == 0) {
          vec.message = paste0("Cell Count = ", 2)
          return(vec.message)
        }
      }
    }
    
    kary.table.dcast = dcast(kary.table.sub[, list(MRN, accession_no, col_name)], MRN + accession_no ~ col_name, fun=length)
    if("chromosomal_count" %in% names(kary.table)) {
      kary.table.dcast[, chromosomal_count := toString(unique(kary.table[, (chromosomal_count)]))]
    }
  }
  if("clone_no" %in% names(kary.table) & is.data.table(kary.table)) {
    kary.table.dcast[, clone.count := length(unique(kary.table[, clone_no]))]
    if(exists("full.deletion.insufficient")){
      kary.table.dcast[, clone.count := (length(unique(kary.table[, clone_no])) - length(full.deletion.insufficient))]
    }
  }
  if("ploidy" %in% names(kary.table) & is.data.table(kary.table)) {
    kary.table.dcast[, ploidy := toString(unique(kary.table[, (ploidy)]))]
  }
  if("cp" %in% names(kary.table) & is.data.table(kary.table)) {
    kary.table.dcast[, cp := toString(unique(kary.table[, (cp)]))]
  }
  if("NA.NA.NA" %in% names(kary.table.dcast) & exists("kary.table.dcast")) {
    kary.table.dcast[,NA.NA.NA := NULL]
    kary.table.dcast[, Normal := T]
    
  }
  return(kary.table.dcast)
}




#' Total Cell Count
#'
#' Returns total cell counts (summed) from pathology text
#'
#' @param text.example Pathology Report Text
#' @param extract.report.text `text.example` is from a Pathology Report (Default: `TRUE`)
#' @return cell.total
#' @examples cell.total = cell_count_total(text.example)
#' @export
#'

cell_count_total <- function(text.example, extract.report.text  = T)  {
  if(extract.report.text) {
    text.example = as.character(text.example)
    text.example = cleanPath.str (path.str = text.example)
    text.example = text.example[grepl("Modal karyotype[:]|Modal karyotype[(]s[)][:]", text.example)]
    if (length(text.example) == 0) {
      return(NA)
      break
    }
    
    if (grepl("Modal karyotype[:]", text.example)) {
      text.example = gsub(".*Modal karyotype[:]", "", text.example)  #Isolate Karyotype, left
    } else {
      text.example = gsub(".*Modal karyotype[(]s[)][:]", "", text.example)  #Isolate Karyotype, left
    }
    stop.words = c("Cytoscan.*|ANALYSIS.*|FISH.*|nuc.*|TEST.*|DISCLAIMER.*|DIAGNOSTIC.*|TestPerformed.*|METAPHASE.*|e-mailed.*|BONE.*|Note.*|Interim.*|PREPARATION.*|Probe.*|yes.*|ABNORMAL.*|not.*|\\*.*|random.*|Int(.)?.*")
    
    #  text.example = gsub("Cytoscan.*|ANALYSIS.*|FISH.*|nuc.*|TEST.*|DISCLAIMER.*|DIAGNOSTIC.*|e-mailed.*|BONE.*|Note.*|Interim.*|PREPARATION.*|Probe.*|yes.*|ABNORMAL.*|not.*|Int(.)?.*", "", text.example)
    text.example = gsub(stop.words, "", text.example)
    
    text.example = gsub('^/', '', text.example) #remove leading blank
    text.example = gsub('\\/', '', text.example) #remove slash
    text.example = gsub("\\n", "", text.example, perl = T) # Remove Carriage Returns
    text.example = gsub("\\\\n", "", text.example)
    text.example = gsub("idemdel", "idem,del", text.example) # Remove Typos
    text.example = gsub("\\s", "", text.example) # Remove White Spaces
    text.example = gsub("\\?", "", text.example) # Remove Question Mark
    
    
    text.example = unlist(strsplit(text.example, "]")) # Split on square bracket , unlist
    text.example = paste0(text.example, "]") # Add back bracket
    text.example = setdiff(text.example, "]")
    text.example = setdiff(text.example, ".]")
    if (length(text.example) == 0) { # For bad reports :(
      return(NA)
    } else {
      
      lst.text = split(text.example, 1:length(text.example)) # split karyotype into list
      
      cell.counts = gsub("\\[|\\]", "", regmatches(lst.text, gregexpr("\\[.*?\\]", lst.text))) # Extract Cell Counts
      cell.counts = gsub("cp", "", cell.counts) # Swap out CP
      cell.counts = as.integer(cell.counts)
      cell.total = sum(cell.counts)
      return(cell.total)
    }
  } else {
    
    text.example = gsub('^/', '', text.example) #remove leading blank
    text.example = gsub('\\/', '', text.example) #remove slash
    text.example = gsub("\\n", "", text.example, perl = T) # Remove Carriage Returns
    text.example = gsub("\\\\n", "", text.example)
    text.example = gsub("idemdel", "idem,del", text.example) # Remove Typos
    text.example = gsub("\\s", "", text.example) # Remove White Spaces
    text.example = gsub("\\?", "", text.example) # Remove Question Mark
    
    
    text.example = unlist(strsplit(text.example, "]")) # Split on square bracket , unlist
    text.example = paste0(text.example, "]") # Add back bracket
    text.example = setdiff(text.example, "]")
    text.example = setdiff(text.example, ".]")
    if (length(text.example) == 0) { # For bad reports :(
      return(NA)
    } else {
      
      lst.text = split(text.example, 1:length(text.example)) # split karyotype into list
      
      cell.counts = gsub("\\[|\\]", "", regmatches(lst.text, gregexpr("\\[.*?\\]", lst.text))) # Extract Cell Counts
      cell.counts = gsub("cp", "", cell.counts) # Swap out CP
      cell.counts = as.integer(cell.counts)
      cell.total = sum(cell.counts)
      return(cell.total)
    }
  }
}


#' @import data.table
NULL

#' Clean Text
#'
#' Removes spaces and symbols from pathology string
#'
#' @param path.str pathology report string
#'
#' @return substring
#'
#' @export
#'


cleanPath.str <- function(path.str) {
  str = sub('^\\s+', '', path.str) #remove leading blank
  str = gsub('[[:blank:]]+', ' ', str) #consolidate blanks
  str = gsub('\r\n', '\n', str) # convert mac formatted newlines to normal newlines
  str = gsub('[[:blank:]]?[\r\n][[:blank:]]?', '\n', str) #remove leading/trailing blanks
  str = gsub('[[:blank:]]?_+[[:blank:]]?', '', str) #remove line breaks of repeating dashes
  str = gsub('[[:blank:]]?\\.\\.+[[:blank:]]?', '', str) #remove line breaks of repeating periods
  str = gsub(':\\s+', ': ', str) # remove extra formatting line breaks after colons
  str = gsub('(ANALYSIS):?[\r\n]+', '\\1: ', str) # make ANALYSIS headers uniform
  str = gsub('\\.[[:blank:]]*Comment', '[\r\n][\r\n]Comment', str) # comments in new paragraph.
  str = gsub('blastd', 'blast', str) # spelling mistake.
  
  return(str)
}


#' Tokenize Text
#'
#' Makes string vector into words
#'
#' @param path.str pathology report string
#'
#' @return vector of words
#'
#' @export
#'

tokenizePath.vec <- function(path.str) {
  pathClean.str = cleanPath.str(path.str)
  vec.path = unlist(strsplit(pathClean.str, '\n\n+'))
  
  return(vec.path)
}

#' @import data.tree
#' @import data.table
NULL

#' Format abnormalities from parsed output
#' @param dt.sample
#'
#' @return mat.sample
#' @export

getFeatureMatrix <- function(dt.sample) {
  dt.parents = dt.sample[,.(clone.ID, parent_clone)]
  clone_IDS = dt.sample$clone.ID
  dt.parents[, parent_clone := gsub("Clone\\s*", "", parent_clone)]
  dt.parents[, parent_clone := as.integer(parent_clone)]
  dt.parents = dt.parents[!is.na(parent_clone)]
  mat.sample = as.matrix(dt.sample[, -c("parent_clone", "Normal")], rownames = "clone.ID")
  if(nrow(dt.parents) > 0) {
    for (i in 1:nrow(dt.parents)) {
      parent.clone.ID = dt.parents[i, parent_clone]
      daughter.clone.ID = dt.parents[i, clone.ID ]
      mat.sample[as.character(daughter.clone.ID), ] = mat.sample[as.character(daughter.clone.ID), ] + mat.sample[as.character(parent.clone.ID), ]
    }
  }
  features = names(colSums(mat.sample)[colSums(mat.sample) > 0])
  if (length(features) > 1) {
    mat.sample = mat.sample[,c(features)]
    rownames(mat.sample) = paste0("Clone_", clone_IDS)
  } else {
    mat.sample = mat.sample[,c(features)]
    mat.sample = as.matrix(mat.sample, ncol = 1)
    rownames(mat.sample) = paste0("Clone_", clone_IDS)
    colnames(mat.sample) = features
  }
  return(mat.sample)
}


#' Function to produce a list of trees
#'
#' @param lst.clones.parsed List of clones parsed into matrices
#' @param vec.accession.no Accession No of interest contained in `lst.clones.parsed`
#'
#' @return lst.unit.trees
#' @export

getCytoTree = function(lst.clones.parsed, vec.accession.no) {
  lst.unit.trees = list()
  for (i in vec.accession.no) {
    print(i)
    if (!i %in% names(lst.clones.parsed)) {
      next
    } else if (nrow(lst.clones.parsed[[i]]) == 1) {
      next
    } else {
      mat.example = getFeatureMatrix(lst.clones.parsed[[i]])
      #      mat.example[mat.example > 1] = 1
      if(ncol(mat.example) > 1) {
        mat.consolidated.full.sub = do.call(cbind, lapply(unname(split.default(as.data.frame(mat.example), apply(mat.example, 2, paste, collapse = ''))),   function(x) matrix(x[,1], dimnames = list(NULL, paste(colnames(x), collapse='\n')))))
        rownames(mat.consolidated.full.sub) = rownames(mat.example)
        mat.consolidated.full.sub[mat.consolidated.full.sub > 1] = 1
        print(mat.consolidated.full.sub)
        mat.consolidated.full = as.matrix(mat.consolidated.full.sub[,names(sort(colSums(mat.consolidated.full.sub), decreasing  = T))])
        colnames(mat.consolidated.full) = names(sort(colSums(mat.consolidated.full.sub), decreasing  = T))
      } else {
        mat.consolidated.full  = mat.example
        mat.consolidated.full[mat.consolidated.full > 1] = 1
        mat.consolidated.full = as.matrix(mat.consolidated.full[,names(sort(colSums(mat.consolidated.full), decreasing  = T))], ncol = 1)
        colnames(mat.consolidated.full) = colnames(mat.example)
      }
      if(ncol(mat.consolidated.full) > 1) {
        clone.complete = names(rowSums(mat.consolidated.full)[rowSums(mat.consolidated.full) >  0])
        mat.consolidated.full = mat.consolidated.full[c(clone.complete),] # Remove empty rows
      }
      dt.example = as.data.table(mat.consolidated.full, keep.rownames = "Clone.ID")
      mlt.example = melt(dt.example, id.vars = 'Clone.ID')
      mlt.example = mlt.example[value > 0]
      mlt.example = mlt.example[order(value, decreasing = T)]
      dt.example[, pathString := mlt.example[, paste(sort(variable), collapse = "/"), by = 'Clone.ID']$V1]
      dt.example[, pathString := paste0("Normal/", pathString)]
      test.tree <- as.Node(data.frame(dt.example[,.(Clone.ID, pathString)]), pathDelimiter = "/")
      lst.unit.trees[[i]] = test.tree
    }
  }
  return(lst.unit.trees)
}


#' Function to get tree dimensions
#'
#' @param tree
#'
#' @return dt.dim
#' @export

getTreeDim = function(tree) {
  df.tree = ToDataFrameTree(tree, "Clone.ID", "leafCount", "level")
  dt.tree = as.data.table(df.tree)
  dt.dim = data.table(height = max(dt.tree$level), width = max(dt.tree$leafCount))
  return(dt.dim)
}

