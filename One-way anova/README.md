## MixTwice on one-way ANOVA

MixTwice (Zheng et al, 2021, Bioinformatics) is a large-scale hypothesis testing tool by variance mixing. It examines the test statistic (t stat for example) as the ratio of effect size and standard error of effect size. As an empirical bayes implementation, MixTwice involves the estimation of those two mixing distributions.

However, MixTwice can only allows the test statistic to be t stat. In other words, the problem setting is only based on two-group comparison testing equal means or one-group test for zero mean.

This follow-up study considers for each test units as a one-way ANOVA (multiple group compairson) where a F stat is always used to compute p-value.