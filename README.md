# Nozzle-exercise
A code for the exercise on the monodimensional and isentropic nozzle simulation

Il codice è stato sviluppato con Matlab 2023a
Richiede la presenza di 3 file:
  La geometria dell’ugello
  I dati per il confronto fra il modello monodimensionale ed il CFD
  I file si chiamano:
      «nozzle_geometry.mat»
      «cfd_data_freeslip.mat»
      «cfd_data_noslip.mat»
Il codice non richiede nessun toolbox aggiuntivo (Matlab base)

Il codice utilizza funzionalità standard per caricare dati (load()), scegliere fra varie condizioni (if…else) e raffigurare I risultati (plot())
L’operazione più «avanzata» che viene effettuata è la risoluzione numerica di una equazione in una variabile (risoluzione numerica di 𝑀=𝑓(𝐴))
Per la risoluzione viene utilizzata fzero()
  https://it.mathworks.com/help/matlab/ref/fzero.html
  https://en.wikipedia.org/wiki/Brent's_method
fzero() richiede di definire una funzione (M_i=f(A_i)) e cerca la soluzione intorno ad un punto x_0 specificato dall’utente o, come nel nostro caso, all’interno di un intervallo specificato, ad es. (0,1) o (1,5)


![immagine](https://github.com/user-attachments/assets/3d94479b-9098-413a-8ae6-c5dc7d98b822)
