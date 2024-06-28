"""
Definition of urls for GVMDF.
"""

from django.urls import path, include
from app import forms, views
from rest_framework import routers
from rest_framework.urlpatterns import format_suffix_patterns

#router = routers.DefaultRouter()
#router.register(r'species', views.species_viewset, basename = 'species')
#router.register(r'protein', views.protein_viewset, basename = 'protein')

urlpatterns = format_suffix_patterns([
     path('', views.api_root),
     path('home/', views.homepage.as_view(), name = 'homepage'),
     path('species/', views.species_list.as_view(), name = 'species_list'),
     path('ls/', views.ls_list.as_view(), name = 'ls_list'),
     path('proteins/', views.protein_list.as_view(), name = 'protein_list'),
     path('country/', views.country_list.as_view({ 'get': 'list' }), name = 'country_list'),
     path('ridges/', views.ridge_plot.as_view(), name = 'ridge_plot'),
     path('classes/', views.mutation_classes.as_view(), name = 'mutation_classes'),
     path('single/', views.single_gene.as_view(), name = 'single_gene'),
     path('global/', views.global_profile.as_view(), name = 'global_profile'),
     path('upset/', views.upset_venn.as_view(), name = 'upset_venn'),
     path('pcr/', views.primer_design.as_view(), name = 'primer_design')
     ])
     #path('species/<int:pk>/', species_get, name = 'species_get'),
     #path('protein/', protein_list, name = 'protein_list'),
     #path('protein/<int:pk>/', protein_detail, name = 'protein_detail')
    

#urlpatterns = [
#    path('', views.home, name='home'),
#    path('contact/', views.contact, name='contact'),
#    path('about/', views.about, name='about'),
#    path('login/',
#         LoginView.as_view
#         (
#             template_name='app/login.html',
#             authentication_form=forms.BootstrapAuthenticationForm,
#             extra_context=
#             {
#                 'title': 'Log in',
#                 'year' : datetime.now().year,
#             }
#         ),
#         name='login'),
#    path('logout/', LogoutView.as_view(next_page='/'), name='logout'),
#    path('admin/', admin.site.urls),
#]
