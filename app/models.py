"""
Definition of models.
"""

from django.db import models
from django.utils.translation import gettext_lazy as _

# Create your models here.

class Country(models.Model):
    class Regions(models.TextChoices):
        EUROPE = "EU", _("Europe")
        ASIA = "AS", _("Asia")
        NORTHAMERICA = "NA", _("North America")
        SOUTHAMERICA = "LA", _("South America")
        OCEANIA = "OC", _("Oceania")
        AFRICA = "AF", _("Africa")

    country = models.CharField(max_length = 32, primary_key = True)
    region = models.CharField(max_length = 2, choices = Regions.choices, db_index = True)

    def __str__(self):
        return self.country + '|' + self.region

class Species(models.Model):
    species = models.TextField()
    refseq = models.TextField(null = True)
    annotation = models.TextField(null = True)
    seglength = models.TextField(null = True)

    def __str__(self):
        return self.species + '|' + self.refseq + '|' + self.annotation

class Lineage(models.Model):
    lineage = models.CharField(max_length = 16, primary_key = True)
    species = models.ForeignKey(Species, db_index = True, on_delete=models.CASCADE, default = 1)

    def __str__(self):
        return self.lineage+ '|' + str(self.species.species)

class Sample (models.Model):
    id = models.CharField(max_length = 32, primary_key = True)
    segment = models.PositiveIntegerField(default = 0)
    sample = models.TextField()
    time = models.DateField(db_index = True, null = True)
    country = models.ForeignKey(Country, db_index = True, on_delete = models.SET_NULL, null = True)
    lineage = models.ForeignKey(Lineage, db_index = True, on_delete = models.SET_NULL, null = True)
    species = models.ForeignKey(Species, db_index = True, on_delete = models.CASCADE, null = True)

    def __str__(self):
        return self.id + '|' + str(self.segment) + '|' + self.sample + '|' + str(self.time) + '|' + self.country.country + '|' + self.lineage.lineage + '|' + self.species.species

class Protein(models.Model):
    protein = models.CharField(max_length = 32, db_index = True)
    annotation = models.TextField()
    species = models.ForeignKey(Species, db_index = True, on_delete=models.CASCADE, default = 1, null = True)
    length = models.PositiveIntegerField(default = 0)
    segment = models.PositiveIntegerField(default = 0)

    def __str__(self):
        return str(self.id) + '|' + self.protein + '|' + self.annotation + '|' + self.species.species + "|" + str(self.length) + "|" + str(self.segment) + "|" + str(self.nested)

class Position(models.Model):
    id = models.PositiveBigIntegerField(primary_key=True, db_index = True)
    rpos = models.PositiveIntegerField()
    rvar = models.TextField()
    qvar = models.TextField()
    qpos = models.PositiveIntegerField()
    sample = models.ForeignKey(Sample, db_index = True, on_delete = models.CASCADE)
    protein = models.ForeignKey(Protein, db_index = True, on_delete = models.SET_NULL, null = True)

    def __str__(self):
        return str(self.id) + '|' + str(self.rpos) + '|' + self.rvar + '|' + self.qvar + '|' + str(self.qpos) + '|' + self.sample.id + '|' + self.protein.protein
    
class Extra(models.Model):
    class Varclass(models.TextChoices):
        NoArg = "NA", _("NoArg")
        SNP = "SN", _("SNP")
        SNP_STOP = "ST", _("SNP Stop")
        SNP_silent = "SS", _("SNP Silent")
        deletion = "DE", _("Deletion")
        deletion_frameshift = "DF", _("Deletion Frameshift")
        deletion_stop = "DT", _("Deletion Stop")
        extragenic = "EX", _("Extragenic")
        insertion = "IN", _("Insertion")
        insertion_frameshift = "IF", _("Insertion Frameshift")
        insertion_stop = "IT", _("Insertion Stop")

    id = models.PositiveBigIntegerField(primary_key=True, db_index = True)
    M_type = models.TextField()
    PM_type = models.TextField()
    pro_variant = models.TextField()
    variant = models.TextField()
    varclass = models.CharField(max_length = 2, choices = Varclass.choices)
    position = models.ForeignKey(Position, db_index = True, on_delete = models.CASCADE, null = True)

    def __str__(self):
        return str(self.id) + '|' + self.M_type + '|' + self.PM_type + '|' + self.pro_variant + '|' + self.variant + '|' + self.varclass + '|' + str(self.position.id)
