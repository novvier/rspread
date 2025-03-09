#' Calcula la ponderación de los niveles de presión sonora
#'
#' Esta función aplica diferentes ponderaciones (A, C o Z) a los niveles de presión sonora
#' en distintas bandas de frecuencia, basándose en los resultados proporcionados.
#'
#' @param rsl Una lista de resultados con niveles de presión sonora en bandas de frecuencia específicas.
#'        Debe contener los elementos `f63`, `f125`, `f250`, `f500`, `f1000`, `f2000`, `f4000` y `f8000`.
#' @param pond Un carácter que indica el tipo de ponderación a aplicar. Puede ser:
#'        - "A": Ponderación A, que simula la respuesta del oído humano a niveles bajos.
#'        - "C": Ponderación C, que tiene una respuesta más plana, utilizada en niveles altos.
#'        - "Z": Ponderación Z (plana), sin modificaciones en la respuesta de frecuencia.
#' @param all_rsl Un valor lógico. Si es `TRUE`, se espera que `rsl` contenga dos listas de resultados
#'        y se usará la lista `noise_propagation`. Si es `FALSE`, se espera que `rsl` contenga ocho listas de resultados.
#'
#' @return Un valor numérico que representa el nivel de presión sonora ponderado.
#'
#' @examples
#' rsl <- list(f63 = 50, f125 = 55, f250 = 60, f500 = 65, f1000 = 70, f2000 = 75, f4000 = 80, f8000 = 85)
#' weighting(rsl, pond = "A")
#'
#' @export
weighting <- function(rsl, pond = "A", all_rsl = F){
  if (all_rsl) {
    if(length(rsl) != 2){
      stop("Se requieren dos listas de resultados para calcular la ponderación")
    }
    rsl <- rsl$noise_propagation
  } else {
    if(length(rsl) != 8){
      stop("Se requiere ocho listas de resultados para calcular la ponderación")
    }
  }
  if (pond == "A") {
    f63 <- 10^((rsl$f63 - 26.2)/10)
    f125 <- 10^((rsl$f125 - 16.1)/10)
    f250 <- 10^((rsl$f250 - 8.6)/10)
    f500 <- 10^((rsl$f500 - 3.2)/10)
    f1000 <- 10^((rsl$f1000 + 0)/10)
    f2000 <- 10^((rsl$f2000 + 1.2)/10)
    f4000 <- 10^((rsl$f4000 + 1.0)/10)
    f8000 <- 10^((rsl$f8000 - 1.1)/10)
    Lp <- 10 * log10(f63 + f125 + f250 + f500 + f1000 + f2000 + f4000 + f8000)
    return(Lp)
  } else if (pond == "C") {
    f63 <- 10^((rsl$f63 - 0.8)/10)
    f125 <- 10^((rsl$f125 - 0.2)/10)
    f250 <- 10^((rsl$f250 + 0)/10)
    f500 <- 10^((rsl$f500 + 0)/10)
    f1000 <- 10^((rsl$f1000 + 0)/10)
    f2000 <- 10^((rsl$f2000 - 0.2)/10)
    f4000 <- 10^((rsl$f4000 - 0.8)/10)
    f8000 <- 10^((rsl$f8000 - 3.0)/10)
    Lp <- 10 * log10(f63 + f125 + f250 + f500 + f1000 + f2000 + f4000 + f8000)
    return(Lp)
  } else if (pond == "Z") {
    f63 <- 10^(x$f63/10)
    f125 <- 10^(rsl$f125/10)
    f250 <- 10^(rsl$f250/10)
    f500 <- 10^(rsl$f500/10)
    f1000 <- 10^(rsl$f1000/10)
    f2000 <- 10^(rsl$f2000/10)
    f4000 <- 10^(rsl$f4000/10)
    f8000 <- 10^(rsl$f8000/10)
    Lp <- 10 * log10(f63 + f125 + f250 + f500 + f1000 + f2000 + f4000 + f8000)
    return(Lp)
  }
}
