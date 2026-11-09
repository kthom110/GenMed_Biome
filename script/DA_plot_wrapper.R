############################################
#############Libraries######################
############################################
library(tidyr)
library(dplyr)
#library(tibble)
library(stringr)
library(phyloseq)
library(ggplot2)
library(ggpmisc)
library(ggrepel)


############################################
#############Prep Function##################
############################################
da_plt_prep = function(obj, type = "auto", 
                       feat.w = "Feature", feat.clean = "Feature", feat.full = "Feature", feat.def = "Definition",
                       coef.col = "coef", pval = "p.raw", signifcol = "qval",
                       cuts = c(0.05, 0.1, 0.2), label.crit = c("x < 0.05"), ntop = Inf,
                       middle = 0){
  
  if(!(obj$input$FeatureType %in% c("Taxonomic", "Functional"))){
    print("Requested Analysis is not DA")
    return(NA)
  }
  
  res.df = obj$output$results
  
  if(type == "auto"){
    toget = c()
    toselect = c("Feature_W", "Feature_Clean", "Feature_Full", "Feature_Definition", "Coefficient", "pval", "padj")
    if(obj$input$DAType == "ZicoSeq"){
      coef.col = obj$input$OOI
      pval = "p.raw"
      signifcol = "qval"
      toget = c(toget, coef.col, pval, signifcol)
      middle = 0
    }
    if(obj$input$DAType == "CoxPH"){
      coef.col = "HR"
      pval = "p.raw"
      signifcol = "qval"
      toget = c(toget, coef.col, pval, signifcol)
      middle = 1
    }
    if(obj$input$FeatureType == "Taxonomic"){
      toget = c(toget, "Feature", "Definition")
      df.tmp = res.df %>% select(all_of(toget)) %>% 
        mutate("Feature_W" = Feature,
               "Feature_Clean" = substr(Feature, start = 4, stop = nchar(Feature)) %>%
                 str_replace_all(pattern = paste0(c("-", "_"), collapse = "|"), replacement = " "),
               "Feature_Full" = Feature,
               "Feature_Definition" = Definition,
               "Coefficient" = c_across(all_of(coef.col)),
               "pval" = c_across(all_of(pval)),
               "padj" = c_across(all_of(signifcol)))
      df.tmp = df.tmp %>% select(all_of(toselect))
    }
    if(obj$input$FeatureType == "Functional"){
      toget = c(toget, "Feature", "Definition", "Original")
      df.tmp = res.df %>% select(all_of(toget)) %>% 
        mutate("Feature_W" = Feature,
               "Feature_Clean" = Feature,
               "Feature_Full" = Original,
               "Feature_Definition" = Definition,
               "Coefficient" = c_across(all_of(coef.col)),
               "pval" = c_across(all_of(pval)),
               "padj" = c_across(all_of(signifcol)))
      df.tmp = df.tmp %>% select(all_of(toselect))
      
      if(obj$input$FeatureSubtype %in% c("EC", "KO")){
        df.tmp = df.tmp %>% mutate(Feature_Clean = Feature_Full)
      } 
      if(obj$input$FeatureSubtype %in% c("Pathway")){
        df.tmp = df.tmp %>% mutate(Feature_Clean = str_split(string = Feature_Full, pattern = ":", simplify = T) %>% 
                                     data.frame() %>% 
                                     select(1) %>% 
                                     unlist(),
                                   Feature_Definition = str_split(string = Feature_Full, pattern = ":", simplify = T) %>% 
                                     data.frame() %>% 
                                     select(2) %>% 
                                     unlist() %>% str_trim())
      } 
    }
  } else{
    toget = c(feat.w, feat.clean, feat.full, feat.def, coef.col, pval, signifcol)
    toget = sapply(toget, function(x){return(which(colnames(res.df) == x))})
    df.tmp = res.df[,toget]
    colnames(df.tmp) = c("Feature_W", "Feature_Clean", "Feature_Full", "Feature_Definition", "Coefficient", "pval", "padj")
  }
  
  #get significance column and direction
  cut = sort(cuts, decreasing = F)
  cut.expr = sapply(1:length(cut), function(i){
    if(i == 1){
      op = paste0("padj < ", cut[i])
    } else{
      lower.tmp = cut[i-1]
      upper.tmp = cut[i]
      op = paste0("padj >= ", lower.tmp, " & padj < ", upper.tmp)
    }
    return(op)
  })
  cut.nm = gsub(pattern = ">=", replacement = "\u2265", x = cut.expr)
  cut.expr = gsub(pattern = "padj", replacement = "df.tmp$padj", x = cut.expr)
  
  
  cut.col = sapply(1:length(cut.expr), function(i){
    exp.tmp = cut.expr[i]
    nm.tmp = cut.nm[i]
    vec.tmp= eval(parse(text = exp.tmp))
    vec.tmp = ifelse(vec.tmp == T, yes = nm.tmp, no = NA)
    return(vec.tmp)
  })
  cut.col = apply(cut.col, MARGIN = 1, function(x){
    op.tmp = x %>% na.omit()
    if(length(op.tmp) == 0){
      op.tmp = "ns"
    }
    return(op.tmp)
  })
  df.tmp = df.tmp %>% mutate(signif= cut.col,
                             direction = case_when(Coefficient >= middle ~ "Positive",
                                                   Coefficient < middle ~ "Negative"))
  df.tmp = df.tmp %>% 
    mutate(signif_level = paste0(direction, " & ", signif)) %>%
    mutate(signif_level = case_when(signif == "ns" ~ "ns", .default = signif_level))
  lvl.vec = c(paste0("Positive & ", cut.nm), paste0("Negative & ", cut.nm), "ns")
  df.tmp = df.tmp %>% 
    mutate(signif_level = factor(signif_level, levels = lvl.vec)) %>%
    select(-c(signif, direction))
  
  #label features that need to be labeled
  if(is.null(label.crit)){
    df.tmp = df.tmp %>% mutate(label_val = "")
  } else{
    label.crit = gsub(pattern = "x", replacement = "df.tmp$padj", x = label.crit)
    lab.tf = ifelse(eval(parse(text=label.crit)), yes = T, no = F)
    df.tmp = df.tmp %>% mutate(label_val = lab.tf)
    
    if(!is.infinite(ntop)){
      df.tmp = df.tmp %>% arrange(pval)
      df.tmp[(ntop+1):nrow(df.tmp),"label_val"] = F
    }
    df.tmp = df.tmp %>% 
      mutate(label_val = case_when(label_val == T ~ Feature_Clean,
                                   .default = ""))
    
  }
  
  #if using "auto", change "Coefficient" column to "HR" if CoxPH
  if(type == "auto"){
    if(obj$input$DAType == "CoxPH"){
      df.tmp = df.tmp %>% dplyr::rename(HR = Coefficient)
    }
  }
  
  #output
  op = list("obj" = obj, "df" = df.tmp)
  return(op)
}

############################################
#############Plot Function##################
############################################
da_plt = function(obj.prep, type = "auto", coef.col = "Coefficient", p.col = "pval", model.type = "ZicoSeq",
                  pt.size = 3, stroke.size=0, pt.color = "auto",  
                  xlab.size = "auto", ylab.size = "auto",
                  xlab.param = "auto", ylab.param = "auto",
                  label.size.param = "auto",
                  tick.size = "auto", text.size='auto',
                  tick.color = "black",
                  label.color = "black", box.pad.size = 0.35, point.pad.size = 0.5,
                  seg.color = "gray50", max.overlaps.param = Inf,
                  max.iter.param = 10000, max.time.param = 0.5,
                  force.param = 1, force.pull.param = 1,
                  legend.tf = T, legend.size.param = "auto", legend.title.param = T,
                  skiplabel = F, pt.lim = c(0,1), pt.brks = c(0.0125, 0.025, 0.05, 0.1), pt.rng = c(1,5)){
  plt = NA
  df.tmp = obj.prep$df
  if(type == "auto"){
    model.type = obj.prep$obj$input$DAType
    if("pval" %in% colnames(df.tmp)){
      p.col = "pval"
    }
    if("Coefficient" %in% colnames(df.tmp)){
      coef.col = "Coefficient"
    } else{
      if("HR" %in% colnames(df.tmp)){
        coef.col = "HR"
      }
    }
  }
  
  #plotting
  if(model.type == "ZicoSeq"){
    plt = ggplot(data = df.tmp, aes(x = .data[[coef.col]], y = -log10(.data[[p.col]]), color = signif_level)) 
  }
  if(model.type == "CoxPH"){
    plt = ggplot(data = df.tmp, aes(x = log2(.data[[coef.col]]), y = -log10(.data[[p.col]]), color = signif_level))
  }
  if(is.numeric(pt.size)){
    plt = plt + geom_point(size = pt.size, stroke = stroke.size) + theme_classic() +
      scale_x_continuous(limits = symmetric_limits)
  } else{
    plt = plt + geom_point(aes(size = .data[[pt.size]]), stroke = stroke.size) + theme_classic() +
      scale_x_continuous(limits = symmetric_limits) + scale_size_continuous(limits = pt.lim, breaks = pt.brks, range = pt.rng)
  }
  
  
  #marker expansion
  #plt = plt + scale_y_continuous(expand = expansion(mult = c(0,0.5)))
  
  #add labels
  if(skiplabel == T | !any(nchar(df.tmp$label_val) > 0)){
    plt = plt
  } else{
    if(label.size.param != "auto"){
      plt = plt +
        geom_text_repel(aes(label = stringr::str_wrap(label_val, 100)), 
                        color = label.color, 
                        box.padding   = box.pad.size,
                        point.padding = point.pad.size,
                        segment.color = seg.color, 
                        max.overlaps = max.overlaps.param,
                        size = label.size.param,
                        min.segment.length = 0,
                        segment.size = unit(0.25, units = "mm"),
                        verbose = T,
                        max.time = max.time.param,
                        max.iter = max.iter.param,
                        force = force.param,
                        force_pull = force.pull.param)
    } else{
      plt = plt +
        geom_text_repel(aes(label = stringr::str_wrap(label_val, 100)), 
                        color = label.color,
                        box.padding   = box.pad.size,
                        point.padding = point.pad.size,
                        segment.color = seg.color, 
                        max.overlaps = max.overlaps.param,
                        min.segment.length = 0,
                        segment.size = unit(0.25, units = "mm"),
                        verbose = T,
                        max.time = max.time.param,
                        max.iter = max.iter.param,
                        force = force.param,
                        force_pull = force.pull.param)
    }
  }
  #alter theme components if necessary
  if(xlab.param != "auto"){
    plt = plt + xlab(label = paste0(xlab.param))
  }
  if(ylab.param != "auto"){
    plt = plt + ylab(label = paste0(ylab.param))
  }
  if(xlab.size != "auto"){
    plt = plt + theme(axis.title.x = element_text(size = xlab.size))
  }
  if(ylab.size != "auto"){
    plt = plt + theme(axis.title.y = element_text(size = ylab.size))
  }
  if(text.size != "auto"){
    plt = plt + theme(axis.text = element_text(size = text.size))
  }
  if(tick.color != "auto"){
    plt = plt + theme(axis.text = element_text(color = tick.color))
  }
  if(tick.size != "auto"){
    plt = plt + theme(axis.text = element_text(size = tick.size))
  }
  if(legend.tf == F){
    plt = plt + theme(legend.position = "none")
  }
  if(legend.title.param == F){
    plt = plt + theme(legend.title = element_blank())
  }
  if(legend.size.param != "auto"){
    plt = plt + theme(legend.text = element_text(size = legend.size.param))
  }
  
  
  if(length(pt.color) != 1){
    plt = plt + scale_color_manual(values = pt.color)
  } else{
    if(pt.color == 1 & pt.color != "auto"){
      plt = plt + scale_color_manual(values = pt.color)
    }
  }
  
  return(plt)
  
}
