#####################################################
# Plot Trinidad geo and rivers...
lib<-c("cowplot","ggplot2","data.table","viridis","patchwork","sf","ggpubr","raster")
lapply(lib,library,character.only=T)
source("~/Exeter/five_aside/five_aside_functions.R")

# Read in the file again
trinidad <- st_read(dsn = "data/trinidad_tobago/trinidad_tobago_coastline_osm_20141020.shp")
# rivers <- st_read(dsn = "data/trinidad_tobago/five_aside_rivers.shp")
all_rivers <- st_read(dsn = "data/trinidad_tobago/all_rivers_trinidad.shp")
aripo_river <- st_read(dsn = "data/trinidad_tobago/aripo_river_trinidad.shp")
guanapo_river <- st_read(dsn = "data/trinidad_tobago/guanapo_river_trinidad.shp")
tacarigua_river <- st_read(dsn = "data/trinidad_tobago/tacarigua_river_trinidad.shp")
madamas_river <- st_read(dsn = "data/trinidad_tobago/madamas_river_trinidad.shp")
oropouche_river <- st_read(dsn = "data/trinidad_tobago/oropouche_river_trinidad.shp")

# # Get elevation data
# trini_elevation <- raster("data/TTO_alt/TTG_alt")
# trini_elevation_df <- as.data.frame(trini_elevation, xy = TRUE)

trini_elevation1 <- raster("data/TTO_alt/n10_w061_1arc_v3.tif")
trini_elevation2 <- raster("data/TTO_alt/n10_w062_1arc_v3.tif")
trini_elevation1_df <- as.data.frame(trini_elevation1, xy = TRUE)
colnames(trini_elevation1_df)[3] <- "alt" 
trini_elevation2_df <- as.data.frame(trini_elevation2, xy = TRUE)
colnames(trini_elevation2_df)[3] <- "alt"
trini_elevation_df <- rbind(trini_elevation1_df,trini_elevation2_df)

# Plot snippet with rivers
trini_elevation_df$x2 <- abs(trini_elevation_df$x)
trini_elevation_df$y2 <- abs(trini_elevation_df$y)

# Full trinidad + elevation...
full <- ggplot() +
  geom_sf(data = trinidad)+
  theme_void()+
  xlim(c(62,60.8))+
  ylim(c(10,10.9))+
  theme(panel.border = element_rect(colour="black",fill=NA,size=3))

# Trim
trini_elevation_df <- trini_elevation_df[trini_elevation_df$x2 > 60.9 &
                                            trini_elevation_df$x2 < 61.7,]
trini_elevation_df <- trini_elevation_df[trini_elevation_df$y2 > 10.55 &
                                           trini_elevation_df$y2 < 10.9,]

# Just elevation...
ggplot() +
  geom_point(data = trini_elevation_df[trini_elevation_df$alt > 0,], aes(x = x, y = y, colour = alt),shape=15)+
  coord_equal(ratio=1/cos(mean(trini_elevation_df$x2)*pi/180))


# plot
river_plot<-ggplot() +
  geom_sf(data = trinidad,colour="white")+
  geom_tile(data = trini_elevation_df[trini_elevation_df$alt > 0,], aes(x = x, y = y, colour = alt,fill=alt,size=0.6),show.legend = F)+
  scale_fill_viridis_c(option = "D")+
  scale_color_viridis_c(option = "D")+
  #scale_colour_gradient(low = "red2", high = "black")+
  #scale_fill_gradient(low = "gray75", high = "black")+
  xlim(c(-61.7,-60.9))+
  ylim(c(10.55,10.9))+
  #coord_equal(ratio=1/cos(mean(trini_elevation_df$x2)*pi/180))
  geom_sf(data = all_rivers,colour="lightblue")+
  geom_sf(data = tacarigua_river,colour=five_aside_colour_rivers$colour[1],size=2)+
  geom_sf(data = guanapo_river,colour=five_aside_colour_rivers$colour[2],size=2)+
  geom_sf(data = aripo_river,colour=five_aside_colour_rivers$colour[3],size=2)+
  geom_sf(data = oropouche_river,colour=five_aside_colour_rivers$colour[4],size=2)+
  geom_sf(data = madamas_river,colour=five_aside_colour_rivers$colour[5],size=2)+
  theme_void()

# Remove mountains...
river_plot2<-ggplot() +
  geom_sf(data = trinidad)+
  #geom_tile(data = trini_elevation_df[trini_elevation_df$alt > 0,], aes(x = x, y = y, colour = alt,fill=alt,size=0.6),show.legend = F)+
  #scale_fill_viridis_c(option = "D")+
  #scale_color_viridis_c(option = "D")+
  #scale_colour_gradient(low = "red2", high = "black")+
  #scale_fill_gradient(low = "gray75", high = "black")+
  xlim(c(-61.7,-60.9))+
  ylim(c(10.55,10.9))+
  #coord_equal(ratio=1/cos(mean(trini_elevation_df$x2)*pi/180))
  geom_sf(data = all_rivers,colour="blue2",alpha=0.5)+
  geom_sf(data = tacarigua_river,colour=five_aside_colour_rivers$colour[1],size=2)+
  geom_sf(data = guanapo_river,colour=five_aside_colour_rivers$colour[2],size=2)+
  geom_sf(data = aripo_river,colour=five_aside_colour_rivers$colour[3],size=2)+
  geom_sf(data = oropouche_river,colour=five_aside_colour_rivers$colour[4],size=2)+
  geom_sf(data = madamas_river,colour=five_aside_colour_rivers$colour[5],size=2)+
  theme_void()

# Make random plot to steal legend
random_plot <- five_aside_colour_rivers
random_plot$meh <- 5:1
random_plot$meh2 <- 1:5

legend_plot <- ggplot(random_plot,aes(meh,meh2,colour=river_F))+
  geom_line(size=2)+
  scale_colour_manual(breaks = five_aside_colour_rivers$river_F,
                      values = five_aside_colour_rivers$colour)+
  # guides(color=guide_legend(override.aes=list(fill=NA)))+
  labs(colour="River")+
  theme(legend.key=element_blank(),
        legend.text = element_text(size=16),
        legend.title = element_text(size=18),
        legend.position = c(0.9,0.8),
        legend.justification = "right")

river_legend <- get_legend(legend_plot)

# Plot together
plot.with.inset <-
  ggdraw() +
  draw_plot(river_plot2) +
  draw_plot(full, x = 0.1, y = .6, width = .3, height = .4)+
  draw_plot(river_legend,x=0.45,y=0.5,width=.3,height=.4)

plot.with.inset_no_mountains <-
  ggdraw() +
  draw_plot(river_plot) +
  draw_plot(full, x = 0.1, y = .6, width = .3, height = .4)+
  draw_plot(river_legend,x=0.5,y=0.5,width=.3,height=.4)

# Save as an R object
saveRDS(plot.with.inset,
        "outputs/five_aside_rivers.rds")

# Save mountains figure...
pdf("figs/FigureSX_river_mountains.pdf",width=15,height=10)
plot.with.inset_no_mountains
dev.off()

