# seidh 1.0
### First release of Seidh library
>En el antiguo idioma nórdico, seiðr (seidh), era un tipo de brujería o hechizo practicado en la edad de hierro escandinava con la finalidad de predecir el futuro. Esta práctica fue asociada a menudo con Odín, por su vasto conocimiento atemporal.

Las interacciones entre proteínas y péptidos representan una gran fracción de todos los procesos bioquímicos que tienen lugar en cualquier organismo vivo. Hoy en día, en pleno desarrollo de la medicina personalizada, el conocimiento sobre cómo se acoplan las proteínas y los péptidos resulta de vital importancia para el desarrollo de nuevos tratamientos y terapias encaminadas a la cura de enfermedades genéticas, oncológicas, autoinmunes y un largo etcétera.

La predicción de estructuras tridimensionales del complejo proteína-péptido (también proteína-proteína) es un sector muy concreto de la biología computacional y estructural. Implica un profundo conocimiento de las interacciones entre las moléculas implicadas, y cómo estas pueden afectar a las estructuras secundaria y terciaria de los complejos; además de un potente banco de recursos bioinformáticos para el manejo de esta información.

La librería Seidh implementa un algoritmo basado en el vuelo de Cuculus canorus (cuco común), un ave parásita que en época de puesta invade nidos de otros pájaros para depositar sus huevos. Este movimiento, que alterna movimientos cortos seguidos de vuelos más largos, resulta particularmente interesante para la predicción de estructuras.

La librería también incluye otros módulos auxiliares para llevar a cabo su cometido. Algunos de estos módulos son de recomendado acceso antes de la invocación de la función principal. Para ello, es recomendable leer detenidamente los requisitos de uso de Seidh y cómo pueden estos otros módulos ayudarte para obtener la estructura del péptido problema.

## Requisitos y guía rápida de uso
Seidh 1.0 está escrito para Python 2.7 y requiere los siguientes módulos de Python:
*	biopython 1.66 o superior
*	numpy 1.14.0 o superior
*	PeptideBuilder 1.0.3 o superior
Para utilizarla, basta con importar el módulo principal:
```{python}
from seidh import seidh
```
Para llamar la función principal solo es necesario indicar los objetos BioPython PDB que contengan la proteína y el binding point, junto con un string que contenga la secuencia FASTA.
```{python}
seidh ( protein, peptide, fasta )
```
Adicionalmente, también puede especificar el número de núcleos de procesador con los que ejecutar la simulación. Por ejemplo, para una simulación que emplee cuatro núcleos:
```{python}
seidh ( protein, peptide, fasta, NUM_PROC=4 )
```
Es preciso notar que Seidh incluye en su librería auxiliar funciones para facilitar la tarea, como open_fasta_peptide(), que independientemente de que se le proporcione un fichero .FAS o una secuencia FASTA, la función comprueba que sea una secuencia válida y devuelve un string conteniendo la secuencia no-ambigua. Para más información, consulta la documentación.
Seidh lleva a cabo una predicción “fine-grained”. En su configuración básica, creará hasta 50 péptidos de prueba, que clasificará y ordenará en función de su calidad. Tenga en cuenta que puedes obtener varios péptidos como resultado, en caso de que Seidh determine que existen varias conformaciones posibles.
