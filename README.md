# Costa-et-al.-2023
Publication: Modularity in host-parasite mixed-networks

## Data: 
1. Input data to construct networks:

Fish.parasite.csv

LOC: river sectors

Fish.sp: fish species sampled in the Guaraguaçu river

COD: sample code

N.IND: number of each individual

Cod.ind: code of each individual

Grupo.parasito: taxon group of each parasite

Gen.parasito: genera of each parasite

Eppt.parasito: eppt of each parasite

Species: parasite species

ABD: parasite abundance in each host individual

LAT: latitude

LONG: longitude

2. Parasite species traits

traits.p.csv

Group: taxon group of each parasite

Family: parasite family

Species: parasite species

ABD: total abundance of each parasite species in the Guaraguaçu river

Dev.stage: development stage found in the Guaraguaçu river of each parasite species

Habitat: aquatic environment in which parasite species is commonly found

Transmition.mode: form of transmission of each parasite species

FinalHost: the final host of each parasite

Life.Cycle: complexity of each parasite life cycle 

Parasitism.form: parasitism form of each parasite species

Organ: which organ the parasite was found in the host

2. Fish individual traits

traits.csv

COD.COL: sample code

N.IND: number of each individual

COD.IND: code of each individual

DATA: day-month-year of the sample collection

SETOR: river sector

RD.CV: gillnet size

ESP: host species

COMP.T.cm: total length

COMP.P.cm: standart length

PESO.g.: host weight

SEXO: host sex

EST.GONADA: gonadal stage, mad: mature; imat: imature; esg: spent phase ; dsv: under development

Parasito: if the host was infected by at least one parasite - sim; if host wasn't infected by at least one parasite - nao

Habitat: host habitat

SW: swimming factor of each host species

SW.CAT: swimming category of each host species

2. Host condiction factor
kn.host.csv

COD.IND: code of each individual

ESP: host species

Kn: individual conditional factor

Kn.l: logarithmized conditional factor

## Code:
Analises_Modularity_form.R

Construction and analyses of mixed networks modularity.
