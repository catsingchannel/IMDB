"""
Definition of views.
"""

from django.http import Http404, FileResponse
from django.shortcuts import render
from rest_framework import status, mixins, generics, viewsets
from rest_framework.decorators import api_view
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework.reverse import reverse
from rest_framework.parsers import JSONParser
from rest_framework.renderers import JSONRenderer
from app.models import *
from app.serializers import *
import app.utils as utils
import os

@api_view(['GET'])
def api_root(request, format = None):
    return Response({
        'home' : reverse('homepage', request = request, format = format),
        'species' : reverse('species_list', request = request, format = format),
        'country' : reverse('country_list', request = request, format = format),
        'lineage and segment' : reverse('ls_list', request = request, format = format),
        'protein' : reverse('protein_list', request = request, format = format),
        'ridge plot' : reverse('ridge_plot', request = request, format = format),
        'mutation classes' : reverse('mutation_classes', request = request, format = format),
        'single gene' : reverse('single_gene', request = request, format = format),
        'global profile' : reverse('global_profile', request = request, format = format)
        })


class homepage(APIView):
    def get(self, request):
        return FileResponse(open(os.path.join(os.getcwd(), 'app', 'templates', 'app', 'homepage.html'), 'rb'))

class species_list(APIView):
    def get(self, request, format = None):
        return Response(SpeciesSerializer(Species.objects.all().values('species'), many = True).data)

class ls_list(APIView):
    def post(self, request, format = None):
        try:
            pk = Species.objects.filter(species = request.data.get('species'))[0]
        except:
            return Response({ 'err' : request.data.get('species') })

        res = { 'segment' : False, 'lineageList' : [], 'seglength' : pk.seglength }

        #if(Sample.objects.filter(species = pk)[0].segment != 0):
        if(pk.id in [3, 4]):
            res["segment"] = True
            res["segmentList"] = [1,2,3,4,5,6,7,8]
        #    for it in Sample.objects.values('segment').filter(species = pk).distinct():
        #        res["segmentList"].append(it["segment"])'

        for it in Lineage.objects.values('lineage').filter(species = pk):
            res["lineageList"].append(it["lineage"])

        return Response(res)

class protein_list(APIView):
    def post(self, request, format = None):
        try:
            pk = Species.objects.filter(species = request.data.get('species'))[0].id
        except:
            return Response({ 'err' : request.data.get('species') })

        return Response(Protein.objects.filter(species__id = pk).values_list('protein', flat = True))

class ridge_plot(APIView):
    def post(self, request, format = None):
        try:
            res = utils.country_ridge_plot(request.data)
        except Exception as e:
            return Response({ 'err' : str(e) })
        return Response(res)

class mutation_classes(APIView):
    def post(self, request, format = None):
        try:
            res = utils.mutation_classes(request.data)
        except Exception as e:
            return Response({ 'err' : str(e) })

        return Response(res)

class single_gene(APIView):
    def post(self, request, format = None):
        try:
            res = utils.mutation_rate(request.data)
        except Exception as e:
            return Response({ 'err' : str(e) })

        return Response(res)

class global_profile(APIView):
    def post(self, request, format = None):
        try:
            res = utils.global_profile(request.data)
        except Exception as e:
            return Response({ 'err' : str(e) })

        return Response(res)

class primer_design(APIView):
    def post(self, request, format= None):
        #try:
        res = utils.assays(request.data)
        #except Exception as e:
        #    return Response({'err':str(e)})

        return Response(res)

class upset_venn(APIView):
    def post(self, request, format = None):
        try:
            res = utils.upset_venn(request.data)
        except Exception as e:
            return Response({'err' : e.__class__.__name__ + str(e)})

        return Response(res)

class country_list(viewsets.ReadOnlyModelViewSet):
    queryset = Country.objects.all()
    serializer_class = CountrySerializer

#class species_viewset(viewsets.ReadOnlyModelViewSet):
#    queryset = Species.objects.all()
#    serializer_class = SpeciesSerializer

#class protein_viewset(viewsets.ReadOnlyModelViewSet):
#    queryset = Protein.objects.all()
#    serializer_class = ProteinSerializer