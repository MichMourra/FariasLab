from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
import urllib.request
import time
import os
import argparse

#############################################################################################
# Para que este programa funcione correctamente es muy importante NO minimizar la ventana del navegador de
# Chrome que se despliegue. Si se requiere realizar otra tarea en ese momento en el equipo, simplemente
# trabar con las otras ventanas enfrente de la del navegador.

__author__ = 'dalopez'

# Parameters:
# 1) --inputPath Path para leer el archivo de entrada (lista).
# 2) --outputPath Path para colocar las carpetas con todas los datos descargados y listas de errores.
# 3) --listNameFile Nombre del archivo con el nombre de todos los datos a descargar.
# 4) --cut Si es igual a "TRUE" se toma unicamente la parte de la secuencia de ENA que corresponde a la
#          secuencia de Scope, de lo contrario (FALSE), se tomara la secuencia completa proporcionada
#          por la misma.

# Ejemplo de ejecucion
# python3 aa_to_nt.py --inputPath . --outputPath . --listNameFile ScopeNewList.txt --cut FALSE

if __name__ == "__main__":
    # Parameter definition
    parser = argparse.ArgumentParser(description='Descarga de secuencias de nts a partir de Scope')
    parser.add_argument("--inputPath", dest="inputPath",
                      help="Path to read input file", metavar="PATH")
    parser.add_argument("--outputPath", dest="outputPath",
                      help="Path to place output folders", metavar="PATH")
    parser.add_argument("--listNameFile", dest="listNameFile",
                      help="List to download", metavar="FILE")
    parser.add_argument("--cut", dest="cut",
                      help="cut the sequence", metavar="NAME",
                      choices=('TRUE', 'FALSE'), default='FALSE')
    args = parser.parse_args()


    # Printing parameter values
    print('-------------------------------- PARAMETROS --------------------------------')
    print("Path para leer el archivo de entrada: " + str(args.inputPath))
    print("Nombre del archivo de lista: " + str(args.listNameFile))
    print("Path de las carpetas finales: " + str(args.outputPath))
    print("cut: " + str(args.cut))


    ##############################################################################################
    START = time.time()

    #Creacion de carpeta contenedora de todas las descargas en caso de que aun no exista.
    newpath = str(args.outputPath) + r'/descargas'
    errorpath = str(args.outputPath) + r'/error'
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    if not os.path.exists(errorpath):
        os.makedirs(errorpath)

    # Obtener toda la lista de cadenas de aa's a buscar
    path_input = str(args.inputPath) + "/" + str(args.listNameFile)
    l = open(path_input, "r")
    lista = l.readlines()

    # Indicar el navegador
    driver = webdriver.Chrome(executable_path="./chromedriver")

    for i in lista:
        print("\n++++++++++++++++++++++++++++++++++++++++++\n")
        i = i.replace("\n", "")
        # Tomar el nombre del registro dado. Ejemplo: "d1tjya1" o "d1jyea_"
        nombre_archivo = i.partition(' ')[0]
        try:

            # Ingresar a la DB SCOPe berkeley
            driver.get("https://scop.berkeley.edu/")
            # Ubicar la caja de texto de busqueda
            secuencia = driver.find_element_by_id("searchbox")
            # Buscar el registro deseado
            secuencia.send_keys(i)
            # Presionar ENTER para realizar la busqueda
            secuencia.send_keys(Keys.ENTER)

            if (str(driver.current_url)[0:32] != 'https://scop.berkeley.edu/sunid='): # En caso de que no sea directo el resultado...
                # Ubicar y dar click en el primer resultado de busqueda
                secuencia = WebDriverWait(driver, 10).until(EC.element_to_be_clickable((By.XPATH, "//*[@id='tab1']/div[1]/div[1]/div/ul[1]/li[1]/a[1]")))
                secuencia.click()

            # Ubicar y dar click en "more details"
            secuencia = WebDriverWait(driver, 10).until(EC.element_to_be_clickable((By.XPATH, "//*[@id='tab1']/div[1]/div[1]/div/div[2]/p[1]/a")))
            secuencia.click()

            if str(args.cut) == 'TRUE':
                # Ubicar cuadro de texto de secuencia y guardarlo en "scope_string"
                scope_string = driver.find_element_by_class_name("seq").text
                # Quitar todos los caracteres inecesarios de la secuencia "scope_string"
                scope_string = scope_string.split('\n')
                scope_string = str(scope_string[1:len(scope_string)])
                scope_string = scope_string.replace('\\n', '')
                scope_string = scope_string.replace('[', '')
                scope_string = scope_string.replace(']', '')
                scope_string = scope_string.replace(',', '')
                scope_string = scope_string.replace(' ', '')
                scope_string = scope_string.replace('\'', '')
                scope_string = scope_string.upper()

            # Ubicar y dar click en el link de Uniprot
            secuencia = driver.find_element_by_xpath("//*[@id='tab1']/div[1]/div[1]/div/ul/li[1]/ul/li/a")
            secuencia.click()

            if str(args.cut) == 'TRUE':
                # Ubicar cuadro de texto de secuencia y guardarlo en "uniprot_string"
                uniprot_string = driver.find_element_by_class_name("sequence").text
                # Quitar todos los caracteres inecesarios de la secuencia "uniprot_string"
                uniprot_string = uniprot_string.split('\n')
                uniprot_string = str(uniprot_string)
                uniprot_string = uniprot_string[9:len(uniprot_string)]
                uniprot_string = uniprot_string.replace('\\n', '')
                uniprot_string = uniprot_string.replace('[', '')
                uniprot_string = uniprot_string.replace(']', '')
                uniprot_string = uniprot_string.replace(',', '')
                uniprot_string = uniprot_string.replace(' ', '')
                uniprot_string = uniprot_string.replace('\'', '')
                for numero in range(10):
                    uniprot_string = uniprot_string.replace(str(numero), '')

                # determinar el inicio y fin de la secuencia de nucleotidos que corresponde unicamente
                # a la secuencia original de scope
                inicio_corte = (uniprot_string.find(scope_string[0:20])*3)
                final_corte = (uniprot_string.find(scope_string[0:20]) + len(scope_string))*3


            #Ubicar y dar click en el segundo link de secuencia de nucleotidos
            ubicado = 0
            try: # En caso de que el boton este en la primer tabla
                secuencia = driver.find_element_by_xpath("//*[@id='cross_references']/table[1]/tbody/tr[1]/td[2]/a[2]")
                driver.execute_script("arguments[0].click();", secuencia)
                ubicado = 1
            except:
                pass

            if ubicado == 0: # En caso de que el boton este en la segunda tabla
                secuencia = driver.find_element_by_xpath("//*[@id='cross_references']/table[2]/tbody/tr[1]/td[2]/a[2]")
                driver.execute_script("arguments[0].click();", secuencia)


            # Ubicar el boton de descargar FASTA y copiar la direccion de enlace
            url_download = WebDriverWait(driver, 10).until(EC.element_to_be_clickable((By.XPATH, "//*[@id='view-container']/div[2]/div/div[2]/div[2]/span[2]/a"))).get_attribute('href')
            # descargar el archivo y guardar en carpeta descargas
            path_guardar = newpath + "/" + '_' + nombre_archivo + ".fasta"
            urllib.request.urlretrieve(url_download, path_guardar)


            # Reabrir archivo recien descargado y eliminar saltos de lineas ('\n')
            with open(path_guardar, "r") as f:
                lineas = f.readlines()
                titulo = lineas[0]

            lineas2 = str(lineas[1:len(lineas)]).replace('\\n', '')
            lineas2 = lineas2.replace('[', '')
            lineas2 = lineas2.replace(']', '')
            lineas2 = lineas2.replace(',', '')
            lineas2 = lineas2.replace(' ', '')
            lineas2 = lineas2.replace('\'', '')


            # Recortar la secuencia de nucleotidos si es el caso
            if str(args.cut) == 'TRUE':
                # Verificar que si se encontro la secuencia de scope dentro de la de uniprot antes de recortar
                if uniprot_string.find(scope_string[0:20]) != -1:

                    lineas2 = lineas2[inicio_corte:final_corte]
                    print("Secuencia CORTADA")

                else:
                    print("Secuencia NO CORTADA")


            # Guargar la secuencia final sin saltos de linea
            with open(newpath + "/" + nombre_archivo + ".fasta", 'w') as f:
                f.write(titulo + lineas2)


            os.remove(path_guardar)
            print("COMPLETE: " + i)


        except:
            cadena = "ERROR: " + i
            print(cadena)
            errorpath_complete = errorpath + '/' + 'error_' + str(args.listNameFile)
            errores = open(errorpath_complete, "a")
            errores.write(i + '\n')
            errores.close()


    driver.close()

        #secuencia.send_keys("aeriafipklvgvgfftsggngaqeagkalgidvtydgptepsvsgqvqlvnnfvnqgydaiivsavspdglcpalkramqrgvkiltwdsdtkpecrsyyinqgtpkqlgsmlvemaahqvdkekakvaffyssptvtdqnqwvkeakakisqehpgweivttqfgyndatkslqtaegiikaypdldaiiapdanalpaaaqaaenlkrnnlaivgfstpnvmrpyvqrgtvkefglwdvvqqgkisvyvanallknmpmnvgdsldipgigkvtvspnseqgyhyeakgngivllpervifnkdnidkydf")
        #time.sleep(3)
        #secuencia = driver.find_element_by_xpath("//*[@id='b1']")
        #secuencia.click()
        #time.sleep(2)

        #1er boton
        #"//*[@id="sequences"]/table[1]/tbody/tr[1]/td[2]/a[1]"
        #"//*[@id="cross_references"]/table[1]/tbody/tr[1]/td[2]/a[1]"
        #"//*[@id="cross_references"]/table[1]/tbody/tr[1]/td[2]/a[1]"
        #2do boton
        #"//*[@id="sequences"]/table[1]/tbody/tr[1]/td[2]/a[2]"
        #"//*[@id="cross_references"]/table[1]/tbody/tr[1]/td[2]/a[2]"
        #"//*[@id="cross_references"]/table[1]/tbody/tr[1]/td[2]/a[2]"
        #3ro boton
        #"//*[@id="sequences"]/table[1]/tbody/tr[1]/td[2]/a[3]"
        #"//*[@id="cross_references"]/table[1]/tbody/tr[1]/td[2]/a[3]"
        #"//*[@id="cross_references"]/table[1]/tbody/tr[1]/td[2]/a[3]"
        #4to boton
        #"//*[@id="sequences"]/table[1]/tbody/tr[1]/td[2]/a[4]"
        #"//*[@id="cross_references"]/table[1]/tbody/tr[1]/td[2]/a[4]"
        #"//*[@id="cross_references"]/table[1]/tbody/tr[1]/td[2]/a[4]"
    END = time.time()
    print("Tiempo transcurrido (min): " + str(round((END - START)/60, 2)))
