# emmean outputs
emm_lsm <- data.frame(emm$emmeans)
## get x positions
emm_lsm$x <- as.numeric(emm_lsm$time_piont)

emm_ctr <- data.frame(emm$contrasts)

# identify significant values
sig <- subset(emm_ctr, p.value < 0.1)
sig <- left_join(sig, emm_lsm, by='time_point')

# identify the positions for each significant value and txt
sig_df <- data.frame()
sig_txt <- data.frame()
for (i in unique(sig$x)){               
        tmp <- subset(sig, x == i)
        tmp_df <- data.frame(
                a=c(unique(tmp$x)-0.25, unique(tmp$x)-0.25
                        , unique(tmp$x) + 0.25, unique(tmp$x) + 0.25
                        ) 
                , b=c(tmp[tmp$Treatment == 'Control', 'response']+3*tmp[tmp$Treatment == 'Control', 'SE.y']
                        , max(tmp$response+4*tmp$SE.y), max(tmp$response+4*tmp$SE.y)
                        , tmp[tmp$Treatment == 'Treated', 'response']+3*tmp[tmp$Treatment == 'Treated', 'SE.y'] 
                        )
        )
        tmp_df$grp <- i

        tmp_txt <- data.frame(
                x = i, y = max(tmp$response+4*tmp$SE.y) + max(tmp$response+4*tmp$SE.y)/50
                , lab = paste0('p = ', round(unique(tmp$p.value), 3))
                )

        sig_df <- rbind(sig_df, tmp_df)
        sig_txt <- rbind(sig_txt, tmp_txt)
}

# plot
(ggplot(emm_lsm, aes(x=time_point, y=response)) 
        + geom_boxplot(aes(fill=Treatment, middle=response, ymin=lower.CL, ymax=upper.CL, lower=(response-SE), upper=(response+SE))
                , stat='identity', position='dodge'
                ) 
        + theme_classic()       
        + scale_x_discrete(labels=c('1', '2', '3', '4', 'some_point'))
        + scale_color_manual(values=c('black', 'black'))
        + scale_fill_manual(values=c('firebrick', 'steelblue'))
        + labs(y='something')
        + geom_line(data=sig_df, aes(x=a, y=b, group=grp))
        + geom_text(data=sig_txt, aes(x=x, y=y, label = lab), family='serif', size=5)
        + theme(
                axis.title.x=element_blank()
                , axis.title.y=element_text(face='bold', size=18)
                , axis.text.x= element_text(face='bold', size=16)
                , axis.text.y=element_text(size=16)
                , legend.title=element_text(size=18, face='bold')
                , legend.text=element_text(size=16)
                )
)
