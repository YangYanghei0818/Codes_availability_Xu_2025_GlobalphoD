
creat_ml_model <- function (y.data, x.data, var, method = "rf", type = c("regression", 
    "classification"), testing.ratio = 0.3, threads = 60, remove.outlier = TRUE, 
    ...) 
{
    if (class(y.data) != "data.frame") 
        stop("The <y.data> must be a dataframe!")
    if (class(x.data) != "data.frame") 
        stop("The <x.data> must be a dataframe!")
    if (!var %in% colnames(y.data)) 
        paste0("Can not find `", var, "` in <y.data>, please check it!") %>% 
            stop()
    type <- ifelse(length(type) == 1, type, type[1])
    if (type == "regression" && !is.numeric(y.data[, var])) {
        paste0("The `", var, "` in <y.data> should be numeric if the <type> is `regression`!") %>% 
            stop()
    }
    if (type == "classification" && !is.factor(y.data[, var])) {
        paste0("The `", var, "` in <y.data> should be factor if the <type> is `classification`!") %>% 
            stop()
    }
    pred.dat.chk <- lapply(X = seq(ncol(x.data)), function(x) is.numeric(x.data[, 
        x])) %>% unlist()
    if (FALSE %in% pred.dat.chk) 
        stop("All data in <x.data> must be numeric!")
    check.length <- length(unique(rownames(y.data) == rownames(x.data)))
    check.logics <- !unique(rownames(y.data) == rownames(x.data))
    if (check.length > 1 | check.logics) 
        stop("Sample ids in <y.data> can not be matched to those in <x.data>!")
    type <- ifelse(length(type) > 1, type[1], type)
    if (!type %in% c("regression", "classification")) 
        stop("The <type> must be one of `regression` and `classification`!")
    dw.raw <- data.frame(target = y.data[, var], x.data)
    if (type == "regression" & remove.outlier) {
        dw <- lapply(X = colnames(dw.raw), function(variable) {
            dat.tmps <- data.frame(row.names = rownames(dw.raw), 
                vals = dw.raw[, which(colnames(dw.raw) == variable)])
            Q <- quantile(dat.tmps$vals, probs = c(0.25, 0.75), 
                na.rm = FALSE)
            iqr <- IQR(dat.tmps$vals)
            up <- Q[2] + 1.5 * iqr %>% as.vector() %>% as.numeric()
            low <- Q[1] - 1.5 * iqr %>% as.vector() %>% as.numeric()
            dat.tmps$vals[which(dat.tmps$vals < low | dat.tmps$vals > 
                up)] <- NA
            colnames(dat.tmps) <- variable
            dat.tmps
        }) %>% do.call("cbind", .) %>% na.omit()
    }
    else {
        dw <- dw.raw
    }
    set.seed(1234)
    i <- base::sample(nrow(dw), testing.ratio * nrow(dw))
    tests.dat <- as.data.frame(dw[i, ])
    train.dat <- as.data.frame(dw[-i, ])
    if (type == "classification") {
        if (length(unique(train.dat$target)) < 2) 
            stop("At least two classifications are required for trainning dataset!")
        if (length(unique(tests.dat$target)) < 2) 
            stop("At least two classifications are required for testing dataset!")
    }
    cl <- parallel::makePSOCKcluster(threads)
    doParallel::registerDoParallel(cl)
    set.seed(1234)
    if (type == "regression") {
        model <- caret::train(target ~ ., data = train.dat, method = method, 
            ...)
    }
    else {
        model <- caret::train(x = train.dat[, 2:ncol(train.dat)], 
            y = train.dat$target, method = method, ...)
    }
    parallel::stopCluster(cl)
    dat <- list(model = model, model.method = method, model.type = type, 
        tests.dat = tests.dat, train.dat = train.dat, var.name = var)
    return(dat)
}



evaluate_ml_model <- function (model.dat, only.show.class = TRUE) 
{
    if (model.dat$model.type == "regression") {
        model <- model.dat$model
        tests.dat <- model.dat$tests.dat
        train.dat <- model.dat$train.dat
        obs.val <- tests.dat$target
        prd.dat <- tests.dat[, which(colnames(tests.dat) != "target")]
        predictions <- stats::predict(model, newdata = prd.dat)
        rmse <- sqrt(mean((obs.val - predictions)^2))
        r <- psych::corr.test(obs.val, predictions, method = "pearson", 
            adjust = "fdr")$r
        p <- psych::corr.test(obs.val, predictions, method = "pearson", 
            adjust = "fdr")$p
        plotdata <- data.frame(obs = obs.val, prd = predictions)
        labedata1 <- data.frame(x = -Inf, y = Inf, label = paste0("Training size: ", 
            nrow(train.dat), " | ", "Testing  size: ", nrow(tests.dat)))
        labedata2 <- data.frame(x = -Inf, y = Inf, label = paste0("RMSE = ", 
            round(rmse, 3), ", Pearson R = ", round(r, 3), ", P = ", 
            round(p, 3)))
        ggplot(plotdata, aes(x = obs, y = prd)) + geom_point(shape = 21, 
            size = 3, fill = "gray60", alpha = 0.7) + geom_smooth(method = "lm", 
            formula = y ~ x) + xlab(paste0("#Observed ", model.dat$var.name)) + 
            ylab(paste0("#Predicted ", model.dat$var.name)) + 
            geom_text(data = labedata1, aes(x = x, y = y, label = label), 
                vjust = 2, hjust = -0.02, size = 5, color = "black") + 
            geom_text(data = labedata2, aes(x = x, y = y, label = label), 
                vjust = 4, hjust = -0.02, size = 5, color = "black") + 
            theme_bw() + labs(title = paste0(model.dat$model.method, 
            " model for ", model.dat$var.name)) + theme(axis.text = element_text(size = 12), 
            axis.title = element_text(size = 14), plot.title = element_text(hjust = 0.5, 
                face = "bold"))
    }
    else if (model.dat$model.type == "classification") {
        model <- model.dat$model
        tests.dat <- model.dat$tests.dat
        train.dat <- model.dat$train.dat
        labedata <- data.frame(x = -Inf, y = Inf, label = paste0("Training size: ", 
            nrow(train.dat), " | ", "Testing  size: ", nrow(tests.dat)))
        predicted <- predict(model, newdata = tests.dat[, which(colnames(tests.dat) != 
            "target")], type = "prob")
        label.df <- lapply(colnames(predicted), function(label) as.character(tests.dat$target)) %>% 
            do.call("cbind", .)
        colnames(label.df) <- colnames(predicted)
        for (i in seq(ncol(label.df))) label.df[, i] <- ifelse(label.df[, 
            i] == colnames(label.df)[i], 1, 0) %>% as.numeric()
        label.df <- as.data.frame(label.df)
        rownames(label.df) <- rownames(tests.dat)
        colnames(predicted) <- paste0(colnames(predicted), "_pred_", 
            model.dat$model.method)
        colnames(label.df) <- paste0(colnames(label.df), "_true")
        final.df <- cbind(label.df, predicted)
        roc.res <- multiROC::multi_roc(final.df, force_diag = F) %>% 
            suppressWarnings()
        plot.roc.df <- multiROC::plot_roc_data(roc.res)
        if (only.show.class) 
            plot.roc.df <- plot.roc.df[which(plot.roc.df$Group %in% 
                unique(tests.dat$target)), ]
        plot.roc.df$Group <- paste0(plot.roc.df$Group, " (AUC = ", 
            round(plot.roc.df$AUC, 3), ")")
        ggplot(plot.roc.df, aes(x = 1 - Specificity, y = Sensitivity)) + 
            geom_path(aes(color = Group), linewidth = 1) + geom_segment(aes(x = 0, 
            y = 0, xend = 1, yend = 1), colour = "grey", linetype = "dotdash") + 
            geom_text(data = labedata, aes(x = x, y = y, label = label), 
                vjust = 2, hjust = -0.02, size = 5, color = "black") + 
            theme_bw() + scale_color_manual(values = RColorBrewer::brewer.pal(9, 
            "Set1")) + labs(title = paste0(model.dat$model.method, 
            " model for ", model.dat$var.name)) + theme(plot.title = element_text(hjust = 0.5, 
            face = "bold"), legend.justification = c(1, 0), legend.position = c(0.999, 
            0.001), legend.title = element_blank(), legend.background = element_rect(fill = NULL, 
            size = 0.5, linetype = "solid", colour = "black"), 
            axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
            legend.text = element_text(size = 12))
    }
}