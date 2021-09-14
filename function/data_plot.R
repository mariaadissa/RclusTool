
### This function plots a given signals dataframe for the given observation (indice) ###

data_plot<-function(signals_data,indice){
fig<-plot_ly()
y_data<-as.data.frame(signals_data[[indice]])
x_data<-c(1:dim(signals_data[[indice]])[1])
var_names<-colnames(y_data)
colnames(y_data)<-var_names

for(i in var_names){ 

fig<- fig %>% add_trace( x=x_data, 
                         y=y_data[,i], 
                         name=i, type="scatter",
                         mode="lines")
}

return(fig)
}



