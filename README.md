# rowcal

rowcal is a r script (soon, a package) for radiocarbon calibration and the analysis of chronological information

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

rowcal operates in the R environment. To use it, you need to install R and the packages 'hdrcde' and 'maptools'


### Installing

To install rowcal, copy the script rowcal.r to your working directory and type 

```
source 'rowcal.r'
```
## Calibrating a date


To calibrate 5000±25 BP

```
cal<-rowcal(5000,25)
```

The object 'cal' can be visualised or summarised thus:

```
plot(cal)
spansc14(cal)
```

### Calculating a KDE for a set of radiocarbon dates

There are many ways of inputing data, the simplest is to copy to the clipboard a two-column list of radiocarbon dates in a spreadsheet. Then:

```
kde<-MCdensity()
plot(kde)
```

## Contributing

email r.mclaughlin@qub.ac.uk

## Versioning

Not yet decided.

## Authors

* **T. Rowan McLaughlin** 

See also the list of [contributors](https://github.com/rowan-mclaughlin/rowcal/contributors) who participated in this project.

## License

This project is licensed under the GNU GENERAL PUBLIC LICENSE Version 3 - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments


