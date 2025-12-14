# TFG: Time Series Fuzzy Clustering Project

Este repositorio contiene la memoria y el código R desarrollado para **clustering difuso de series temporales**, implementando diversos métodos de agrupamiento basados en diferentes métricas de similitud.

## 📋 Descripción

El proyecto se centra en el **análisis cluster de series temporales**, específicamente en la aplicación de técnicas de **clustering difuso** para agrupar series estacionarias generadas por modelos estocásticos similares. El código implementa múltiples enfoques para calcular características y disimilitudes entre series temporales, incluyendo:

1. **Implementación de clustering difuso FCMdC** (Fuzzy C-Medoids Clustering)
2. **Múltiples métricas de distancia** para comparación de series temporales, que podemos agrupar en:
   - **Métricas que discriminan en base a la forma**:
     - **EUCL**: Distancia Euclidiana estándar
     - **DTW**: Dynamic Time Warping
   
   - **Métricas que discriminan en base a la estructura**:
     - **ACF**: Autocorrelaciones
     - **PACF**: Autocorrelaciones parciales
     - **Piccolo**: Coeficientes autoregresivos
     - **QAF**: Autocovarianzas cuantil
      
4. **Simulación de escenarios** con diferentes modelos estocásticos (AR, MA, bilineales, no lineales)
5. **Evaluación de resultados** mediante índices de validación (ARI, Jaccard)

## 📁 Estructura del código

`simulaciones.R` incluye:

- **Funciones auxiliares**: Cálculo de ACF, PACF, Piccolo, QAF
- **Algoritmo FCMdC-QAF**: Implementación principal de clustering difuso
- **Generadores de escenarios**: 7 escenarios diferentes con series temporales simuladas
- **Pipeline completo**: Cálculo de matrices de características, clustering y evaluación

