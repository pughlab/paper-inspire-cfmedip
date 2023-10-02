median_range <- function(x, embellish='') {
  median_string <- value_to_string(median(x, na.rm=T))
  if (grepl('\\.', median_string)) {
    post_decimal_numbers = nchar(strsplit(median_string, '\\.')[[1]][2])
  } else {
    post_decimal_numbers = 0
  }
  
  return(
    sprintf(
      '%s%s [range: %s, %s]',
      median(x, na.rm=T) %>% value_to_string,
      embellish,
      min(x, na.rm=T) %>% round(post_decimal_numbers) %>% sprintf(paste0('%.', as.character(post_decimal_numbers), 'f'), .),
      max(x, na.rm=T) %>% round(post_decimal_numbers) %>% sprintf(paste0('%.', as.character(post_decimal_numbers), 'f'), .)
    )
  )
}

pval_to_string <- function(p, n_digits = 3) {
  threshold = 10^(-n_digits)
  if (p < threshold) {
    return(sprintf('p < %s', threshold))
  } else {
    #return(sprintf('p=%#.3g', p))
    return(sprintf('p=%s', round(p, n_digits)))
  }
}

value_to_string <- function(r) {
  return(sprintf('%#.3g', r))
}

ratio_to_string <- function(r) {
  return(sprintf("%.2f", round(r, 2)))
}

value_to_percent <- function(r) {
  return(paste0(as.character(signif(r, 3) * 100), '%'))
}

get_immune_split_probes <- function(data) {
  immune_nonimmune_split <- data %>%
    distinct(probeID) %>%
    inner_join(data_meth_atlas) %>%
    group_by(probeID) %>%
    summarise(
      immune_mean = mean(beta[immune_cell]),
      nonimmune_mean = mean(beta[!immune_cell]),
      diff_t_value = t.test(beta ~ immune_cell)$statistic,
      diff_p_value = t.test(beta ~ immune_cell)$p.value
    ) %>%
    ungroup() %>%
    arrange(-diff_t_value)
  
  immune_nonimmune_full_included <- data %>%
    distinct(probeID) %>%
    inner_join(immune_nonimmune_split) %>%
    mutate(
      included = ((immune_mean < 0.1) | (immune_mean > 0.6)) & diff_p_value < 0.01,
      immune_group = case_when(
        immune_mean > nonimmune_mean & included ~ 'Immune',
        immune_mean < nonimmune_mean & included ~ 'Non-immune',
        TRUE ~ 'Neither'
      )
    ) %>%
    arrange(immune_group, immune_mean - nonimmune_mean) %>%
    mutate(probeID = factor(probeID, levels = unique(probeID)))
  
  return(immune_nonimmune_full_included)
}

get_ahr_string <- function(coxmodel, coefficient) {
  summary_string = sprintf(
    'aHR = %s (%s)\np = %s',
    signif(exp(coef(coxmodel)[coefficient]), 3),
    signif(exp(confint(coxmodel)[coefficient, ]), 3) %>% paste(collapse = ' - '),
    signif(summary(coxmodel)$coef[coefficient, 'Pr(>|z|)'], 3)
  )
}

get_ahr <- function(coxmodel, coefficient) {
  return(list(
    aHR = signif(exp(coef(coxmodel)[coefficient]), 3),
    lCI = signif(exp(confint(coxmodel)[coefficient, 1]), 3),
    uCI = signif(exp(confint(coxmodel)[coefficient, 2]), 3),
    p = summary(coxmodel)$coef[coefficient, 'Pr(>|z|)'], 3
  ))
}

get_survplot <- function(data, title = '', legend = 'top', legend_title = NULL, cowplot=TRUE, rel_widths = c(2,1), ahrstring = TRUE, ahr_coef_name = 'changeDecrease', color_vector = NULL, getfit = FALSE, ylab = 'Survival Probability', risk.table=T) {
  survmodel <- survfit(Surv(time, event) ~ change, data = data)
  coxmodel <- coxph(Surv(time, event) ~ change + cohort, data = data)
  
  if (is.null(legend_title)) {
    legend_title <- 'Strata'
  }
  
  if (typeof(ahrstring) == 'character') {
    summary_string = ahrstring
  } else if (ahrstring) {
    summary_string = get_ahr_string(coxmodel, ahr_coef_name)
  } else {
    summary_string = ''
  }
  
  names(survmodel$strata) <- gsub("change=", "", names(survmodel$strata))
  
  if (is.null(color_vector)) {
    color_vector = RColorBrewer::brewer.pal(data$change %>% unique %>% length, 'Set1')
  }
  
  survplot <- ggsurvplot(
    survmodel,
    risk.table = risk.table,
    tables.theme = theme_nothing(),
    pval = summary_string,
    pval.coord = c(16, 0.92),
    pval.size = 3.5,
    tables.col = "strata",
    tables.y.text = F,
    break.time.by = 6,
    data = data,
    legend = legend,
    title = title,
    palette = color_vector,
    ylab = ylab
  )
  
  survplot$plot <- survplot$plot + labs(color = legend_title)
  
  survplot_grid <- plot_grid(
    survplot$plot + theme(legend.position = 'none'),
    survplot$table, 
    ncol = 1,
    align = 'v',
    rel_heights = c(80, 20)
  )
  
  if (legend != 'none') {
    survplot_grid <- survplot_grid %>%
      plot_grid(
        get_legend(survplot$plot + labs(color = legend_title)),
        nrow = 1,
        rel_widths = rel_widths
      )
  }
  
  if (cowplot) {
    return_survplot <- survplot_grid
  } else {
    return_survplot <- survplot
  }
  
  if (getfit) {
    output = list(
      survplot = return_survplot,
      fit = survmodel,
      data = data
    )
    if (ahrstring) {
      output$cox = get_ahr(coxmodel, ahr_coef_name)
    }
    return(output)
    
  } else {
    return(return_survplot)
  }
}

get_cohort_specific_survplot <- function(cohort_name) {
  cohort_data <- plotdata_cfmedip_response %>%
    select(
      patient,
      time = pfs_since_C3,
      event = progressed,
      change = change,
      cohort = cohort
    )
  
  if (! cohort_name == 'All patients') {
    cohort_data <- cohort_data %>% filter(cohort == cohort_name)
  }
  
  fit <- survfit(
    Surv(time, event) ~ change,
    data = cohort_data
  )
 
  survplot <- ggsurvplot(
    fit,
    risk.table = TRUE,
    tables.theme = theme_nothing(),
    tables.col = "strata",
    tables.y.text = F,
    break.time.by = 6,
    data = cohort_data,
    title = cohort_name 
  )
  
  survplot_grid <- plot_grid(
    survplot$plot,
    survplot$table, 
    ncol = 1,
    align = 'v',
    rel_heights = c(80, 20)
  )
  
  return(survplot_grid)
}



plot_fragment_diffgram <- function(insert_size_table, comparator_size_table, facet_groups=NULL, color_table=NULL, bootstrap_sample_size = 2000, plot_labels = NULL) {
  if (!is.null(facet_groups)) {
    insert_size_table$cohort = facet_groups
  } else {
    insert_size_table$cohort = 'Test group'
  }
  
  if(is.null(plot_labels)) {
    plot_labels = labs(
      x = 'Insert Size',
      y = 'Proportion Difference',
      color='Group'
    )
  }
  
  plotdata_insert_diff <- insert_size_table %>%
    distinct(cohort, library) %>%
    group_by(cohort) %>%
    summarise(
      cancer_library = sample(library, bootstrap_sample_size, replace = TRUE)
    ) %>%
    mutate(
      normal_library = sample(comparator_size_table$library %>% unique %>% sample(bootstrap_sample_size, replace=TRUE))
    ) %>%
    ungroup() %>%
    left_join(
      insert_size_table %>% select(library, cohort, insert_size, cancer_proportion=proportion), by = c('cancer_library' = 'library', 'cohort') 
    ) %>%
    left_join(
      comparator_size_table %>% select(library, insert_size, normal_proportion=proportion), by = c('normal_library' = 'library', 'insert_size')
    ) %>%
    mutate(
      proportion_diff = cancer_proportion - normal_proportion,
      pair = paste(cancer_library, normal_library)
    ) %>%
    left_join(
      color_table %>% rename(cancer_library = library)
    )
  
  if (is.null(color_table)) {
    plotdata_insert_diff_medians <- plotdata_insert_diff %>%
      group_by(cohort, insert_size) %>%
      summarise(
        median_diff = median(proportion_diff, na.rm=T)
      ) %>%
      ungroup()
    
    line_element <- geom_line(alpha = 0.1, size=0.5)
    groups_overlay <- geom_line(aes(
      y = median_diff,
      group = 1
    ),
    data = plotdata_insert_diff_medians,
    size = 1,
    color = 'black'
    )
  } else {
    plotdata_insert_diff_medians <- plotdata_insert_diff %>%
      group_by(cohort, color_group, insert_size) %>%
      summarise(
        median_diff = median(proportion_diff, na.rm=T)
      ) %>%
      ungroup()
    
    line_element = geom_line(aes(color = color_group), alpha = 0.1)
    
    groups_overlay <- geom_line(aes(
      y = median_diff,
      group = color_group,
      color = color_group
    ),
    data = plotdata_insert_diff_medians,
    size = 1
    )
  }
  
  output_plot <- plotdata_insert_diff %>%
    ggplot(aes(
      x = insert_size,
      y = proportion_diff,
      group = pair
    )) +
    line_element +
    facet_wrap(~cohort) +
    scale_color_viridis_d() +
    groups_overlay +
    xlim(30, 250) +
    plot_labels
  
  return(output_plot)
}