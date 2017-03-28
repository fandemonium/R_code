fit.the.model = function(data, alpha) {

   cpc = vector(length(traits), mode = "list")
   names(cpc) = traits

   # find the parents of each trait (may be genes or other traits).
   for (t in seq_along(traits)) {

     # BLUP away the family structure.
     m = lmer(as.formula(paste(traits[t], "~ (1|FUNNEL:PLANT)")), data = data)
     data[!is.na(data[, traits[t]]), traits[t]] =
       data[, traits[t]] - ranef(m)[[1]][paste(data$FUNNEL, data$PLANT, sep = ":"), 1]
     # find out the parents.
     cpc[[t]] = learn.nbr(data[, c(traits, genes)], node = traits[t], debug = FALSE,
                   method = "si.hiton.pc", test = "cor", alpha = alpha)

   }#FOR

   # merge the relevant variables to use for learning.
   nodes = unique(c(traits, unlist(cpc)))
   # yield has no children, and genes cannot depend on traits.
   blacklist = tiers2blacklist(list(nodes[!(nodes %in% traits)],
                 traits[traits != "YLD"], "YLD"))

   # build the Bayesian network.
   bn = hc(data[, nodes], blacklist = blacklist)

   return(bn)

 }#FIT.THE.MODEL
