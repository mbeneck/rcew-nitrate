# Reynolds Nitrate Project Notes

## Daily Working Notes

### 7-16-26. - AGU abstract data availability

#### New data availability
- RME only has new data for precip currently available to 11/20/2025
- 125 precip data also through 2025 for Dobson


#### Things I Learned
- RME missing samples not as bad as I thought- covers december storms, just not early jan storm
- Master sample sheet got really messed up. had to re-enter formulas for calculating day from julian day for sample time and autosampler tub time. Recalculated and seems to work, but also restored version from 5-18-26. Need to compare to make sure I didn't overwrite any dates


## Todos:
- add script running order to README. I think its scan processing -> Model FIts, Export, Validation?
- fit dobson PLSR 
- import dobson weather, snow, soil temp, soil moisture, discharge
- joining dobson results
turn this data processing into actual functions so its not so manual
- compare old Master Sample Sheet to new, make sure no samples have different dates.