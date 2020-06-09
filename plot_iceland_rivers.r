library('maps')
library('ggplot2')

iceland <- map_data("world", region = "Iceland")

map <- ggplot(iceland, aes(long, lat))+ 
  labs(x=NULL, y=NULL)+
  geom_polygon(fill='grey')+
  geom_path()+ #contour
  theme_bw()+ 
  geom_point(data=riv, aes(long, lat), colour = "red")
  
#map + geom_text(data = riv, aes(long, lat, label = river),
            #size=4, nudge_x =0.05, nudge_y = 0.05,
            #check_overlap=TRUE)
#geom_label removes points

#text overlap so label png
                 