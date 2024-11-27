# rowcal

rowcal is a R package for the analysis of chronological information in archaeology, history and related subjects, including tools for radiocarbon calibration and data analysis. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

rowcal operates in the R environment, there are no additional dependencies. 


### Installing

To install rowcal using the `devtools` package, use: 

```
devtools::install.github('rowan-mclaughlin/rowcal')
```
## Calibrating a date


To calibrate 5000±25 BP

```
cal<-rowcal(5000,25)
```

The object 'cal' can be visualised or summarised thus:

```
plot(cal)
hdr(cal) # 95% highest density region
```

### Calculating a KDE for a set of radiocarbon dates

There are many ways of inputing data, the simplest is to copy to the clipboard a two-column list of radiocarbon dates in a spreadsheet. Then:

```
kde<-MCdensity()
plot(kde)
```

## Contributing

email rowan.mclaughlin@mu.ie

## Versioning

Versions 0.x.x were the orginal R script used in a number of publications since 2016. 
Package first release: 1.0.9000

## Authors

* **T. Rowan McLaughlin** 

See also the list of [contributors](https://github.com/rowan-mclaughlin/rowcal/contributors) who participated in this project.

## License

This project is licensed under the GNU GENERAL PUBLIC LICENSE Version 3 - see the [LICENSE.md](LICENSE.md) file for details

## References
KDE models:
McLaughlin, T.R. 2019. On applications of space-time modelling with open-source 14C age calibration. Journal of Archaeological Method and Theory 26, 479–501. https://doi.org/10.1007/s10816-018-9381-3

Growth rate models:
McLaughlin, T. R., Gómez-Puche, M., Cascalheira, J., Bicho, N. and Fernández-López de Pablo, J. 2021. Late Glacial and Early Holocene human demographic responses to climatic and environmental change in Atlantic Iberia. Philosophical Transactions of the Royal Society B 376, 20190724. https://dx.doi.org/10.1098/rstb.2019.0724

SPDboot:
Fernández-López de Pablo, J., Gutiérres-Roig, M.,Gómez-Puche, M., Silva, F., McLaughlin, R., and Lozano, S. 2019. Palaeo-demographic modelling supports a population bottleneck during the Pleistocene-Holocene transition in Iberia. Nature Communications 10, 1872. http://doi.org/10.1038/s41467-019-09833-3

Calibration curves:
Reimer et al. 2020 \url{doi.org/doi/10.1017/RDC.2020.41}
Heaton et al. 2020 \url{doi.org/doi/10.1017/RDC.2020.68}

MCmoving:
(in prep)

Published literature using these tools:


Brownlee E. (2023) A radiocarbon-based model of changing burial rites in early medieval England. Radiocarbon. 2023;65(6):1232-1252. doi:10.1017/RDC.2023.110

Chapple, R.M., McLaughlin, T.R., and Warren, G. 2022. “…where they pass their unenterprising existence…”: change over time in the Mesolithic of Ireland as shown in radiocarbon dated activity. Proceedings of the Royal Irish Academy: Archaeology, Culture, History, Literature 122C, 1–38. https://doi.org/10.1353/ria.0.0002

Dowd, M., Stimpson, C., Connolly, R., Bonsall, J., Kahlert, T., and McLaughlin, R. 2024. New insights on the fauna of Ireland’s Younger Dryas and Early Holocene from Alice & Gwendoline Cave. Quaternary Science Reviews. https://doi.org/10.1016/j.quascirev.2024.108827

Feeser, I., Dörfler, W., Kneisel, J., Hinz, M., & Dreibrodt, S. (2019). Human impact and population dynamics in the Neolithic and Bronze Age: Multi-proxy evidence from north-western Central Europe. The Holocene, 29(10), 1596-1606. https://doi.org/10.1177/0959683619857223

French, C., Hunt., C.O., Grima, R., McLaughlin., R., Stoddart, S., and Malone, C. (eds.) 2020. Temple Landscapes: Fragility, change and resilience of Holocene environments in the Maltese Islands. McDonald Institute, Cambridge. ISBN 978-1-902937-99-1. https://doi.org/10.17863/CAM.59611 (569 pages).

Gleeson, P. (2022). Reframing the first millennium AD in Ireland: archaeology, history, landscape. Proceedings of the Royal Irish Academy: Archaeology, Culture, History, Literature 122, 87-122. https://dx.doi.org/10.1353/ria.2022.0013.

Gleeson, P., and McLaughlin, R. 2021. Ways of Death: cremation and belief in first-millennium AD Ireland. Antiquity 95, 382–399. https://doi.org/10.15184/aqy.2020.251

Hannah, E. (2021). 800 Years of Making: Early Medieval Craftworking and an 8th-Century Decline in Ireland. Medieval Archaeology, 65(2), 238–268. https://doi.org/10.1080/00766097.2021.1999070

Hannah, E. (2023). A chronology for unenclosed settlements in early medieval Ireland: Settlement patterns in the late first millennium AD. Proceedings of the Royal Irish Academy: Archaeology, Culture, History, Literature 123C. https://doi.org/10.1353/ria.2023.a913615

Hannah, E., and McLaughlin R. 2019. Long-term archaeological perspectives on new genomic and environmental evidence from early medieval Ireland. Journal of Archaeological Science 106, 23–28.  http://doi.org/10.1016/j.jas.2019.04.001 (communicating author)

Hill, E., Lubell, D. and McLaughlin, T. R. 2019. A re-evaluation of radiocarbon dates using Helix melanostoma shell at two Capsian sites – Aïn Misteheyia and Kef Zoura D – In the Télidgène Basin of eastern Algeria. Journal of Archaeological Science: Reports 28, 102004. https://doi.org/10.1016/j.jasrep.2019.102004

Leplongeon, A. (2021). The Main Nile Valley at the End of the Pleistocene (28–15 ka): Dispersal Corridor or Environmental Refugium?. Front. Earth Sci. 8. https://doi.org/10.3389/feart.2020.607183

Mallory, J. 2023. 9 From The Steppe To Ireland: The Impact of aDNA Research. The Indo-European Puzzle Revisited: Integrating Archaeology, Genetics, And Linguistics, P.129.

Malone, C., Grima, R., McLaughlin, R., Parkinson, E. W., Stoddart, S., and Vella, N. (eds). 2020. Temple places: Excavating cultural sustainability in prehistoric Malta. McDonald Institute, Cambridge. ISBN 978-1-913344-03-0. https://doi.org/10.17863/CAM.62630 (851 pages).

McLaughlin, R. 2023. Dating and Chronology. In Hartwell, B., Gormley, S., Brogan, C., and Malone, C. (eds) Ballynahatty: Excavations in a Neolithic Monumental Landscape. Oxford: Oxbow, 149–154. ISBN 9781789259711

McLaughlin, R., Hannah, E., and Coyle-McClung, L. 2018. Frequency analyses of historical and archaeological datasets reveal the same pattern of declining sociocultural activity in 9th to 10th Century CE Ireland. Cliodynamics 9: 1–24. https://doi.org/10.21237/C7clio9136654 

McLaughlin, T. R. 2020. An archaeology of Ireland for the Information Age: A reply to Baillie, Cassidy, Plunkett and Waddell. Emania 25, 61–66.

McLaughlin, T. R. 2020. An archaeology of Ireland for the Information Age. Emania 25, 7–30.

McLaughlin, T. R., Gómez-Puche, M., Cascalheira, J., Bicho, N. and Fernández-López de Pablo, J. 2021. Late Glacial and Early Holocene human demographic responses to climatic and environmental change in Atlantic Iberia. Philosophical Transactions of the Royal Society B 376, 20190724. https://dx.doi.org/10.1098/rstb.2019.0724 

McLaughlin, T.R. 2019. On applications of space-time modelling with open-source 14C age calibration. Journal of Archaeological Method and Theory 26, 479¬¬¬–501. 

McLaughlin, T.R., Whitehouse, N.J., Schulting, R.J., McClatchie, M., Barratt, P., and Bogaard, A. 2016. The changing face of Neolithic and Bronze Age Ireland: a big data approach to the settlement and burial archives. Journal of World Prehistory 29, 117-153.  https://doi.org/10.1007/s10963-016-9093-0

Parkinson, E. W. and McLaughlin, T. R. 2020. Lifeways at the acme of the south Italian Neolithic: new chronological and bioarchaeological data from Masseria Fonteviva, Apulia. Journal of Archaeological Science: Reports 34, 102589. https://doi.org/10.1016/j.jasrep.2020.102589

Parkinson, E. W., McLaughlin, T. R., Esposito, C., Stoddart, S., and Malone, C. 2021. Radiocarbon dated trends and the prehistory of the Central Mediterranean. Journal of World Prehistory 34 (3), 317¬–379. (joint first and communicating author). https://doi.org/10.1007/s10963-021-09158-4 

Parkinson, E.W., Stoddart, S., Sparacello, V. et al. Multiproxy bioarchaeological data reveals interplay between growth, diet and population dynamics across the transition to farming in the central Mediterranean. Sci Rep 13, 21965 (2023). https://doi.org/10.1038/s41598-023-49406-5

Stoddart, S., Power, R., Thompson, J., Mercieca-Spiteri, B., McLaughlin, R., Pace, A. Parkinson, E., and Malone, C. (eds). 2022. Temple People: Bioarchaeology, resilience and culture in prehistoric Malta. McDonald Institute, Cambridge. ISBN 978-1-913344-08-5. https://doi.org/10.17863/CAM.91914 (476 pages).

Thompson, J. E., Parkinson, E. W., McLaughlin, T. R., Barratt, R. P., Power, R. K., Mercieca-Spiteri, B., Stoddart, S., and Malone, C. 2020. Placing and remembering the dead in late Neolithic Malta: bioarchaeological and spatial analysis of the Xagħra Circle Hypogeum, Gozo. World Archaeology 52. https://doi.org/10.1080/00438243.2019.1745680


