#!/usr/bin/env Rscript
#title: "plot_Saknovich"
#author: "Alejandro Lopez"
#date: "7/10/2020"
# Colocar el path del directorio con los outputs arrojados por los programas de Saknovich y asi poder obtener todos 
# los archivos "_elongation_profile.dat"
direccion <- "/home/cmourra/cg_cotrans/output_input_ALL/"
setwd(direccion)
# Obtener la lista de todos los documentos con la terminacion requerida
lista <- list.files(direccion, pattern="_elongation_profile.dat", all.files=FALSE, full.names=FALSE)
if (length(lista) == 0) {
  print("No se encontraron documentos en el directorio indicado")
}else{
  print("Documentos encontrados: ")
  print(lista)
  dir.create("Graps_Saknovich") # Creacion del directorio
}

# Crear y guardar las graficas en imagenes jpg.
for (i in lista) {
  
  tabla <- read.table(i, sep = "\t", header = T)
  
  # tomar la primer parte del nombre del archivo para nombrar asi a la imagen
  name <- strsplit(i, "_")[[1]][1]
  
  # unir con la direccion del directorio de "Graps_Saknovich" y el nombre de la grafica
  direccion_con_nombre <- paste(direccion,"Graps_Saknovich", sep = "/")
  direccion_con_nombre <- paste(direccion_con_nombre,name, sep = "/")
  direccion_con_nombre <- paste(direccion_con_nombre,".png", sep = '')
  
  # Guardar la grafica
  png(file=direccion_con_nombre, width=600, height=400)
  plot(tabla$L,tabla$F_L, xlab = "L (residue)", ylab = "F_L", type= "o", col = "blue", main = i)
  dev.off()
  
  print("Guardado en:")
  print(direccion_con_nombre)
}



