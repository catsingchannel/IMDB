from app.models import *
from rest_framework import serializers

class SpeciesSerializer (serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Species
        fields = ['id', 'species', 'refseq', 'annotation', 'seglength']

class CountrySerializer (serializers.ModelSerializer):
    class Meta:
        model = Country
        fields = ['country', 'region']

#class ProteinSerializer (serializers.HyperlinkedModelSerializer):
#    species = serializers.HyperlinkedRelatedField(view_name = 'species-detail', read_only = True)

#    class Meta:
#        model = Protein
#        fields = ['id', 'protein', 'annotation', 'species']