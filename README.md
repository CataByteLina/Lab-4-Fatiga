# Lab-4-Fatiga
La electromiografía es un estudio que mide la actividad eléctrica de los músculos y los nervios que los controlan. Se usa para diagnosticar trastornos neuromusculares, como neuropatías, miopatías o enfermedades de la motoneurona. El electromiograma es la grabación de dicha actividad electrica. La electromiografía se puede realizar de dos formas, una se realiza insertando electrodos de aguja en los músculos para registrar su actividad en reposo y durante la contracción, la otra se hace superficialmente por medio de electrodos, utilizando 2 electrodos de trabajo y uno de tierra. En este laboratorio se hizo el estudio de el musculo braquiradial de forma superficial.
La respuesta impulsiva es la respuesta del ,usculo o fibra muscular al estimulo, es decir que describe cómo se propaga el potencial de acción en el músculo.

# Procedimiento
1.Se usó el sistema de adquisición de datos DAQ NI USB 6001/6002/6003
2. Se conectaron los electrodos al amplificador y al DAQ, para después ser adheridos en el musculo braquioradial del antebrazo izquierdo sujeto.
3. Se define la frecuencia de muestreo con 1000Hz
4. El sujeto realiza contracciones musculares haciendose un registros de la señal EMG durante este proceso. Una vez capturada la señal se realiza el filtrado.
5. Se realiza el aventanamiento

# Filtro

Se utilizó un pasa banda Butterworth de orden 4 para eliminar el ruido de baja frecuencia como artefactos o línea base y eliminar ruido de alta frecuencia, también porque se centra en el rango de frecuencias normal de las señales electromiograficas.

```def bandpass_filter(signal, fs, lowcut=20, highcut=450, order=4):
    nyq = 0.5 * fs  # Frecuencia de Nyquist
    low, high = lowcut / nyq, highcut / nyq  # Normalización de frecuencias
    b, a = butter(order, [low, high], btype='band')  # Diseño del filtro Butterworth
    return filtfilt(b, a, signal)  # Filtrado con corrección de fase
```

# Aventanamiento

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

Se utilizó la ventana Hamminng porque reduce el efecto de fuga espectral al calcular la Transformada de Fourier y suaviza los bordes de las ventanas para evitar artefactos.  

A su vez se realiza un analisis espectral utilizando la Transformada de Fourier para convertir la señal al dominio de la fecuencia para saber como se distribuyen las frecuencias en la señal capturada. Teoricamente si la frecuencia mediana disminuye

```def spectral_analysis(windows, fs):
    freq_medians = []
    for window in windows:
        fft_vals = np.abs(fft(window))[:len(window)//2]  # Calcula la FFT y toma la mitad positiva
        fft_freqs = np.fft.fftfreq(len(window), 1/fs)[:len(window)//2]  # Calcula las frecuencias correspondientes
        median_freq = fft_freqs[np.where(np.cumsum(fft_vals) >= np.sum(fft_vals)/2)[0][0]]  # Frecuencia mediana
        freq_medians.append(median_freq)
    return freq_medians
```

