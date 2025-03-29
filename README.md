# Lab-4-Fatiga
La electromiografía es un estudio que mide la actividad eléctrica de los músculos y los nervios que los controlan. Se usa para diagnosticar trastornos neuromusculares, como neuropatías, miopatías o enfermedades de la motoneurona. El electromiograma es la grabación de dicha actividad electrica. La electromiografía se puede realizar de dos formas, una se realiza insertando electrodos de aguja en los músculos para registrar su actividad en reposo y durante la contracción, la otra se hace superficialmente por medio de electrodos, utilizando 2 electrodos de trabajo y uno de tierra. En este laboratorio se hizo el estudio de el músculo braquioradial de forma superficial.
La respuesta impulsiva es la respuesta del múusculo o fibra muscular al estimulo, es decir que describe cómo se propaga el potencial de acción en el músculo. 

## Procedimiento

1.Se usó el sistema de adquisición de datos DAQ NI USB 6001/6002/6003

2. Se conectaron los electrodos al amplificador y al DAQ, para después ser adheridos en el músculo braquioradial del antebrazo izquierdo sujeto.
   
3. Se define la frecuencia de muestreo con 1000Hz
 
4. El sujeto realiza contracciones musculares haciendose un registros de la señal EMG durante este proceso. Una vez capturada la señal se realiza el filtrado.
 
5. Se realiza el aventanamiento

6. Análisis espectral

## Filtro

Se aplicó un filtro pasa banda en el rango de 20 a 450 Hz es fundamental para mejorar la calidad de la señal y garantizar la precisión en su procesamiento. El filtro pasa altas, con un umbral de 20 Hz, permite eliminar artefactos de movimiento, fluctuaciones de corriente continua y ruido fisiológico de baja frecuencia, mientras que el filtro pasa bajas, con un límite superior de 450 Hz, atenúa el ruido eléctrico, la interferencia por contacto del electrodo y componentes de alta frecuencia no relevantes. La implementación de este filtrado mediante la función `filtfilt()` asegura un procesamiento sin distorsión de fase, lo que es esencial para la detección precisa de contracciones musculares y el análisis de la fatiga muscular.

![image](https://github.com/user-attachments/assets/610dc76b-f9a4-4989-8c64-0f001fb19502)


```def bandpass_filter(signal, fs, lowcut=20, highcut=450, order=4):
    nyq = 0.5 * fs  # Frecuencia de Nyquist
    low, high = lowcut / nyq, highcut / nyq  # Normalización de frecuencias
    b, a = butter(order, [low, high], btype='band')  # Diseño del filtro Butterworth
    return filtfilt(b, a, signal)  # Filtrado con corrección de fase
```
![image](https://github.com/user-attachments/assets/298151af-cadf-49d9-afd7-c492a073ffc2)

## Aventanamiento

El aventanamiento es una técnica en el procesamiento de señales donde una señal larga se divide en segmentos más pequeños (ventanas). A cada ventana se le aplica una función matemática que suaviza los bordes, evitando problemas como la fuga espectral cuando se hace análisis de frecuencia. Estas ventanas permiten analizar cambios "locales", en este caso es una ventana por cada contracción muscular.

```
def apply_hamming_windows(signal, peaks, fs, window_size=0.2):
    window_samples = int(fs * window_size)  # Convierte el tamaño de la ventana de segundos a muestras
    windows = []
    for peak in peaks:
        start, end = max(0, peak - window_samples // 2), min(len(signal), peak + window_samples // 2)
        if end - start == window_samples:
            windows.append(signal[start:end] * np.hamming(window_samples))  # Se aplica la ventana de Hamming
    return windows
```

Se utilizó la ventana Hamming porque reduce el efecto de fuga espectral al calcular la Transformada de Fourier y suaviza los bordes de las ventanas para evitar artefactos.  

A su vez se realiza un análisis espectral utilizando la Transformada de Fourier para convertir la señal al dominio de la fecuencia para saber como se distribuyen las frecuencias en la señal capturada. Teoricamente si la frecuencia mediana disminuye

```def spectral_analysis(windows, fs):
    freq_medians = []
    for window in windows:
        fft_vals = np.abs(fft(window))[:len(window)//2]  # Calcula la FFT y toma la mitad positiva
        fft_freqs = np.fft.fftfreq(len(window), 1/fs)[:len(window)//2]  # Calcula las frecuencias correspondientes
        median_freq = fft_freqs[np.where(np.cumsum(fft_vals) >= np.sum(fft_vals)/2)[0][0]]  # Frecuencia mediana
        freq_medians.append(median_freq)
    return freq_medians
```

# Análisis Espectral

Se evalua la disminución de la frecuencia de la mediana, si esta disminuye es porque hay fatiga.
```
plt.figure()
plt.plot(freq_medians, marker='o', linestyle='-')
plt.xlabel('Número de ventana')
plt.ylabel('Frecuencia mediana (Hz)')
plt.title('Evolución de la frecuencia mediana')
plt.show()
```
Utilizando la prueba de hipótesis se evalua si el cambio es significativo para poder evidenciar si hubo fatiga o no

```python
def evaluate_fatigue_first_last(freq_medians):
    if len(freq_medians) < 2:
        print(" No hay suficientes contracciones para la prueba de hipótesis.")
        return None, None, "Sin datos", (None, None)
    
    first, last = np.array(freq_medians[0]), np.array(freq_medians[-1])  # Se comparan la primera y la última frecuencia mediana
    
    if len(freq_medians) >= 3:
        _, p_shapiro_first = shapiro(freq_medians[:3])  # Prueba de normalidad en las primeras ventanas
        _, p_shapiro_last = shapiro(freq_medians[-3:])  # Prueba de normalidad en las últimas ventanas
    else:
        p_shapiro_first, p_shapiro_last = None, None
    
    if p_shapiro_first is not None and p_shapiro_first > 0.05 and p_shapiro_last > 0.05:
        t_stat, p_value = ttest_rel([first], [last])  # Prueba T pareada si los datos son normales
        test_name = "Prueba T pareada"
    else:
        t_stat, p_value = wilcoxon([first], [last])  # Prueba de Wilcoxon si los datos no son normales
        test_name = "Prueba de Wilcoxon"
    
    return t_stat, p_value, test_name, (p_shapiro_first, p_shapiro_last)
```
La función`evaluate_fatigue_first_last(freq_medians)` evalúa la fatiga muscular comparando la primera y la última frecuencia mediana de una serie de mediciones de EMG. Primero, verifica si hay suficientes datos; si no, devuelve un mensaje de error. Luego, extrae la primera y la última frecuencia mediana y realiza una prueba de normalidad (Shapiro-Wilk) en las tres primeras y tres últimas mediciones, si hay al menos tres datos disponibles. Dependiendo del resultado de esta prueba, usa una prueba T pareada si los datos son normales o una prueba de Wilcoxon si no lo son. Finalmente, devuelve el estadístico de la prueba, el valor p, el tipo de prueba realizada y los valores p de las pruebas de normalidad. Sin embargo, el código tiene errores porque pasa listas de un solo elemento a las pruebas estadísticas, lo que puede generar fallos. Una mejor versión promedia varias ventanas iniciales y finales antes de aplicar la prueba, asegurando un análisis más robusto y confiable de la fatiga muscular.

![image](https://github.com/user-attachments/assets/c6785897-2ae4-49ea-b1ac-149577514340)

Dado que se comparan la primera y la última contracción, el resultado indica que la mediana de la frecuencia no cambió de manera significativa entre ambas. Esto sugiere que, a lo largo del experimento, no hubo una disminución notable en la frecuencia mediana, lo que generalmente se asocia con la fatiga muscular, esto significa que el músculo mantuvo su comportamiento electromiográfico sin una reducción evidente en la frecuencia, por lo que no se puede concluir que haya ocurrido fatiga en el período analizado.
Esto se puede evidenciar en la grafica de las ventanas.

![image](https://github.com/user-attachments/assets/7b7c31e7-9eb9-4ea5-bf86-723e77802cdd)
![image](https://github.com/user-attachments/assets/5ce8d829-593d-4c75-8762-fbe8f31a5308)

