# rspread

## Descripción
`rspread` es un paquete de R para modelar la propagación sonora de fuentes de área (polígonos) y fuentes puntuales (puntos). Está basado en `SoundMapTools` y permite realizar predicciones precisas utilizando diversas herramientas geoespaciales.

## Instalación
Actualmente, `rspread` no está disponible en CRAN, pero puedes instalarlo directamente desde GitHub con:

```r
# Instalar devtools si no está instalado
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

# Instalar rspread desde GitHub
devtools::install_github("novvier/rspread")
```

## Dependencias
`rspread` requiere R (>= 4.3) y las siguientes dependencias:

- `terra` (>= 1.8)
- `sf` (>= 1.0)
- `dplyr` (>= 1.1)
- `foreign` (>= 0.8)
- `raster` (>= 3.6)
- `gdistance` (>= 1.6.4)

Además, para ejecutar algunas funciones se requiere un compilador de Fortran (ej. `gfortran`).

## Uso básico
Una vez instalado, puedes cargar el paquete y comenzar a usarlo:

```r
library(rspread)
```

Ejemplo de cómo utilizar el paquete para modelar la propagación sonora:

```r
# Ejemplo de código de uso de rspread (por desarrollar)
```

## Reporte de errores
Si encuentras algún error o tienes sugerencias, puedes reportarlas en:
[Repositorio de issues](https://github.com/novvier/rspread/issues)

## Licencia
Este paquete está licenciado bajo la licencia **GPL (>= 2)**.

