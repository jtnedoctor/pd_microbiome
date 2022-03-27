#install.packages("maps")
library(maps)
library(ggplot2)
data("world.cities")
filepath="~project/16s/ndd/allmetadata/merge_metadata.tsv"
worldmetadata <- fread(filepath,sep = '\t',fill = T,quote = F,check.names = F,stringsAsFactors = F)
worldmetadata <- worldmetadata[worldmetadata$diseasetype=="PD"|worldmetadata$diseasetype=="PD_Con"]
design <- worldmetadata
design$project <- toupper(design$project)
design[!duplicated(design$project),] %>%
  select(country,project,platform,region,city) %>%
  dplyr::rename(name=city) -> nddmap
nddcities <- subset(world.cities, name%in%nddmap$name & country.etc%in%nddmap$country)
nddcities <- inner_join(nddcities,nddmap,by='name')
#p <- qplot(long, lat, data = nddcities,colour=project,shape=region,size=1)+ borders("world", size= 0.5)
p <- ggplot()+ borders("world", colour="white", fill="peachpuff", size= 0.5,alpha = 0.5) + #peru,darkseagreen2,peachpuff
  geom_point(data=nddcities,size=2,position=position_jitter(h=0.5, w=0.5),aes(long, lat,color=project,shape=region)) +  #position = "jitter",
  theme(panel.background = element_rect(fill = 'white', colour = 'white'), # theme(panel.background = element_blank())
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) 
p
ggsave(p, file = "ndd_pd_all_rnaseq_16s_worldmap_peachpuff.pdf", width = 17, height = 17, unit = 'cm',bg = "transparent")
