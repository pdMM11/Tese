from django.shortcuts import render
# from django.core.paginator import Paginator
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.pagination import PageNumberPagination
from django_filters.rest_framework import DjangoFilterBackend
from rest_framework import filters
from rest_framework.views import APIView
from rest_framework import generics

from .models import *
from .serializers import *
from django.shortcuts import render
from django.http import Http404, JsonResponse, HttpResponse
import requests
from subprocess import run, PIPE
import sys, os
from urllib.parse import unquote
import pickle
from propythia.machine_learning import MachineLearning
import pandas as pd
from rest_framework import status
from rest_framework.response import Response
from crmapp.ml_models.Propythia_ML_PMML import MachineLearning_PMML
from sklearn2pmml import sklearn2pmml
from pypmml import Model
from datetime import date, datetime
from Bio.Align.Applications import ClustalOmegaCommandline

from weblogo import *


class CustomPagination(PageNumberPagination):
    page_size = 0
    page_size_query_param = 'page_size'
    max_page_size = 50000

    def get_paginated_response(self, data):
        return Response({
            'links': {
                'next': self.get_next_link(),
                'previous': self.get_previous_link()
            },
            'count': self.page.paginator.count,
            'page_size': self.page_size,
            'results': data
        })


def weblogo(request):
    # performs clustal alignment in order to use it in the weblogo analysis
    seqs = unquote(request.GET.get('seq'))
    in_file = "unaligned.fasta"
    file = open("unaligned.fasta", "w")
    file.write(seqs)
    file.close()
    out_file = "out_filename.fasta"
    clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=False)
    print(clustalomega_cline)
    os.system('cmd /c crmapp\clustal-omega-1.2.2-win64\\' + str(clustalomega_cline) + ' --force')

    """
    out_file = "out_filename.clustal_num"
    clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=False)
    print(clustalomega_cline)
    os.system('cmd /c crmapp\clustal-omega-1.2.2-win64\\' + str(clustalomega_cline) + ' --outfmt clustal --force')
    """

    file_out = open("out_filename.fasta", "r")
    seqs_aligned = file_out.readlines()
    # return_data = {'data': seqs_aligned}

    seqs = read_seq_data(file_out)
    logodata = LogoData.from_seqs(seqs)
    logooptions = LogoOptions()
    logooptions.title = "VFP WEBSERVER"
    logoformat = LogoFormat(logodata, logooptions)
    weblogo_txt = txt_formatter(logodata, logoformat)

    weblogo_jpeg = jpeg_formatter(logodata, logoformat)

    weblogo_file = "weblogo.txt"
    weblogo = open(weblogo_file, "w")
    data_weblogo = str(weblogo_txt)[2:len(str(weblogo_txt)) - 1].replace('\\n', '\n').replace('\\t', '\t')
    weblogo.write(data_weblogo)
    weblogo.close()

    file_out.close()
    os.remove(in_file)
    os.remove(out_file)
    # print(return_data)
    output = seqs_aligned[0]
    seq_found = False
    for i in range(1, len(seqs_aligned) - 1):
        if seqs_aligned[i + 1][0] == '>':
            output += seqs_aligned[i]
            seq_found = True
        elif seq_found:
            output += seqs_aligned[i]
            seq_found = False
        else:
            output += seqs_aligned[i][0:len(seqs_aligned[i]) - 1]
    output += seqs_aligned[len(seqs_aligned) - 1]
    return JsonResponse({'data': output}, safe=False)


###########################################################

# Create your views here.

def iedb(request):
    data = {
        'method': unquote(request.GET.get('method')),
        'sequence_text': unquote(request.GET.get('sequence_text')),
        'window_size': unquote(request.GET.get('window_size'))
    }

    response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/bcell/', data=data)
    # return render(request, 'home.html', {'data': response.text})
    # return_data = response.text.split('\n')

    return JsonResponse(response.text, safe=False)


def clustal(request):  # error no urllib

    email = unquote(request.GET.get('email'))
    seqs = unquote(request.GET.get('seq'))

    file = open("unaligned.fasta", "w")
    file.write(seqs)
    file.close()
    # os.system('cmd /c"python crmapp/clustalo.py --email' + email + ' --stype protein' ' unaligned.fasta"')

    command = [sys.executable, 'crmapp/clustalo.py',
               "--email", email,
               "--stype", "protein",
               "unaligned.fasta"]

    output = run(command, shell=True, stdout=PIPE)
    # os.system(command)

    print(output)

    # file = open("unaligned.fasta", "w")
    # file.write(seqs)
    # file.close()
    # os.system('cmd /c"python clustalo.py --email' +  email +  ' --stype protein' ' unaligned.fasta"')

    return render(request, 'home.html')


def ml_predict(request):
    seq = unquote(request.GET.get('sequence'))
    try:
        window_size = int(unquote(request.GET.get('window_size')))
    except:
        window_size = 15
    try:
        gap = int(unquote(request.GET.get('gap')))
    except:
        gap = 1
    try:
        model_picked = unquote(request.GET.get('model'))
    except:
        model_picked = 'svm'

    model = None

    if model_picked == 'svm':
        model = pickle.load(open('crmapp/ml_models/dataset3_all_svm.sav', 'rb'))
        # model = Model.fromFile('crmapp/ml_models/Model_pmml_svm_28052020.pmml')
    elif model_picked == 'rf':
        model = pickle.load(open('crmapp/ml_models/dataset3_all_rf_model.sav', 'rb'))
    elif model_picked == 'gboosting':
        model = pickle.load(open('crmapp/ml_models/dataset3_all_gboosting_model.sav', 'rb'))
    elif model_picked == 'knn':
        model = pickle.load(open('crmapp/ml_models/dataset3_all_knn_model.sav', 'rb'))
    elif model_picked == 'lr':
        model = pickle.load(open('crmapp/ml_models/dataset3_all_lr_model.sav', 'rb'))
    elif model_picked == 'gnb':
        model = pickle.load(open('crmapp/ml_models/dataset3_all_gnb_model.sav', 'rb'))
    elif model_picked == 'nn':
        model = pickle.load(open('crmapp/ml_models/dataset3_all_nn_model.sav', 'rb'))

    dataset = pd.read_csv(r'crmapp/ml_models/dataset3_all_svc.csv', delimiter=',')

    x_original = dataset.loc[:, dataset.columns != 'labels']
    labels = dataset['labels']

    # create Machine learning object
    ml = MachineLearning(x_original, labels, classes=['non_vfp', 'vfp'])
    # ml = MachineLearning_PMML(x_original, labels, classes=['non_vfp', 'vfp'])

    result = None
    if model is not None:
        result = ml.predict_window(model, seq=seq, x=None, window_size=window_size,
                                   gap=gap, features=[], names=None, y=None,
                                   filename=None)

        return JsonResponse(result.to_json(orient='table'), safe=False)

    else:
        JsonResponse([], safe=False)

    """
    if model_picked != 'svm':
        filename = 'crmapp/ml_models/dataset3_all_svc'
        dataset_in = r'crmapp/ml_models/dataset3_all_svc.csv'
        dataset = pd.read_csv(dataset_in, delimiter=',')
        x_original = dataset.loc[:, dataset.columns != 'labels']
        labels = dataset['labels']
        # create Machine learning object
        ml = MachineLearning_PMML(x_original, labels, classes=['non_vfp', 'vfp'])

        best_svm_model = ml.train_best_model('svm')

        # save the model to disk
        # filename2 = filename + 'svm_model.sav'
        # pickle.dump(best_svm_model, open(filename2, 'wb'))
        sklearn2pmml(best_svm_model, "Model_pmml_svm_27052020.pmml", with_repr=True)

        return HttpResponse('<h1>Page was found</h1>')

        return render(request, 'home.html')
        """


class WriteResultsAPIView(APIView):
    @api_view(["POST"])
    @csrf_exempt
    def write_ml_results(request):
        if request.method == "POST":
            data = request.data
            today = datetime.now()
            today = today.strftime("%d-%m-%Y-%H%M%S")
            username = os.getlogin()  # Fetch username
            file = open(f'C:\\Users\\{username}\\Desktop\\results_ml' + str(today) + '.txt', 'w')
            file.write(data['data'])
            file.close()
            return JsonResponse({'response': "Data successfully saved."}, safe=False)
        raise Http404

    @api_view(["POST"])
    @csrf_exempt
    def write_fusion_peptide_results(request):
        if request.method == "POST":
            data = request.data
            today = datetime.now()
            today = today.strftime("%d-%m-%Y-%H%M%S")
            username = os.getlogin()  # Fetch username
            file = open(f'C:\\Users\\{username}\\Desktop\\fusion_peptide_' + str(today) + '.csv', 'w', newline='')
            file.write(data['data'])
            file.close()
            return JsonResponse({'response': "Data successfully saved."}, safe=False)
        raise Http404

    @api_view(["POST"])
    @csrf_exempt
    def write_inhibitor_antibody_results(request):
        if request.method == "POST":
            data = request.data
            today = datetime.now()
            today = today.strftime("%d-%m-%Y-%H%M%S")
            username = os.getlogin()  # Fetch username
            file = open(f'C:\\Users\\{username}\\Desktop\\inhibitor_antibody_' + str(today) + '.csv', 'w', newline='')
            file.write(data['data'])
            file.close()
            return JsonResponse({'response': "Data successfully saved."}, safe=False)
        raise Http404

    @api_view(["POST"])
    @csrf_exempt
    def write_peptide_references_results(request):
        if request.method == "POST":
            data = request.data
            today = datetime.now()
            today = today.strftime("%d-%m-%Y-%H%M%S")
            username = os.getlogin()  # Fetch username
            file = open(f'C:\\Users\\{username}\\Desktop\\peptide_references_' + str(today) + '.csv', 'w', newline='')
            file.write(data['data'])
            file.close()
            return JsonResponse({'response': "Data successfully saved."}, safe=False)
        raise Http404

    @api_view(["POST"])
    @csrf_exempt
    def write_peptide_references_results(request):
        if request.method == "POST":
            data = request.data
            today = datetime.now()
            today = today.strftime("%d-%m-%Y-%H%M%S")
            username = os.getlogin()  # Fetch username
            file = open(f'C:\\Users\\{username}\\Desktop\\peptide_references_' + str(today) + '.csv', 'w', newline='')
            file.write(data['data'])
            file.close()
            return JsonResponse({'response': "Data successfully saved."}, safe=False)
        raise Http404

    @api_view(["POST"])
    @csrf_exempt
    def write_peptide_structures_results(request):
        if request.method == "POST":
            data = request.data
            today = datetime.now()
            today = today.strftime("%d-%m-%Y-%H%M%S")
            username = os.getlogin()  # Fetch username
            file = open(f'C:\\Users\\{username}\\Desktop\\peptide_structures_' + str(today) + '.csv', 'w', newline='')
            file.write(data['data'])
            file.close()
            return JsonResponse({'response': "Data successfully saved."}, safe=False)
        raise Http404

    @api_view(["POST"])
    @csrf_exempt
    def write_protein_results(request):
        if request.method == "POST":
            data = request.data
            today = datetime.now()
            today = today.strftime("%d-%m-%Y-%H%M%S")
            username = os.getlogin()  # Fetch username
            file = open(f'C:\\Users\\{username}\\Desktop\\protein_' + str(today) + '.csv', 'w', newline='')
            file.write(data['data'])
            file.close()
            return JsonResponse({'response': "Data successfully saved."}, safe=False)
        raise Http404

    @api_view(["POST"])
    @csrf_exempt
    def write_protein_references_results(request):
        if request.method == "POST":
            data = request.data
            today = datetime.now()
            today = today.strftime("%d-%m-%Y-%H%M%S")
            username = os.getlogin()  # Fetch username
            file = open(f'C:\\Users\\{username}\\Desktop\\protein_references_' + str(today) + '.csv', 'w', newline='')
            file.write(data['data'])
            file.close()
            return JsonResponse({'response': "Data successfully saved."}, safe=False)
        raise Http404

    @api_view(["POST"])
    @csrf_exempt
    def write_tax_host_results(request):
        if request.method == "POST":
            data = request.data
            today = datetime.now()
            today = today.strftime("%d-%m-%Y-%H%M%S")
            username = os.getlogin()  # Fetch username
            file = open(f'C:\\Users\\{username}\\Desktop\\tax_host_' + str(today) + '.csv', 'w', newline='')
            file.write(data['data'])
            file.close()
            return JsonResponse({'response': "Data successfully saved."}, safe=False)
        raise Http404

    @api_view(["POST"])
    @csrf_exempt
    def write_taxonomy_virus_results(request):
        if request.method == "POST":
            data = request.data
            today = datetime.now()
            today = today.strftime("%d-%m-%Y-%H%M%S")
            username = os.getlogin()  # Fetch username
            file = open(f'C:\\Users\\{username}\\Desktop\\taxonomy_virus_' + str(today) + '.csv', 'w', newline='')
            file.write(data['data'])
            file.close()
            return JsonResponse({'response': "Data successfully saved."}, safe=False)
        raise Http404
        """
        writer = csv.writer(file)
        writer.writerow(['idtaxonomy', 'commonname', 'family', 'genre',
                         'species', 'subspecies', 'ncbitax'])
        for obj in list(queryset.values()):
            writer.writerow(
                [obj['idtaxonomy'], obj['commonname'], obj['family'], obj['genre'],
                 obj['species'], obj['subspecies'], obj['ncbitax']])
        """


class FusionPeptidesAPIView(generics.ListCreateAPIView):
    queryset = FusionPeptides.objects.all()
    serializer_class = FusionPeptidesSerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['idfusion_peptides', 'residues', 'sequence', 'annotation_method',
                        'exp_evidence', 'protein', 'protein__name', 'protein__idtaxonomy',
                        'protein__idtaxonomy__family',
                        'protein__idtaxonomy__commonname']
    search_fields = ['idfusion_peptides', 'residues', 'sequence', 'annotation_method',
                     'exp_evidence', 'protein__name', 'protein__idtaxonomy__family',
                     'protein__idtaxonomy__commonname']

    def get_object(self, pk=None):
        if pk is None: return
        try:
            return FusionPeptides.objects.get(pk=pk)
        except FusionPeptides.DoesNotExist:
            raise Http404

    def put(self, request, pk=None, format=None):
        if pk is None: return Response(status=status.HTTP_400_BAD_REQUEST)
        snippet = self.get_object(pk)
        serializer = FusionPeptidesSerializer(snippet, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk=None, format=None):
        if pk is None: return
        snippet = self.get_object(pk)
        snippet.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class FusionPeptidesAPIView_Save(generics.ListCreateAPIView):
    queryset = FusionPeptides.objects.all()
    serializer_class = FusionPeptidesSerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['idfusion_peptides', 'residues', 'sequence', 'annotation_method',
                        'exp_evidence', 'protein', 'protein__name', 'protein__idtaxonomy',
                        'protein__idtaxonomy__family',
                        'protein__idtaxonomy__commonname']
    search_fields = ['idfusion_peptides', 'residues', 'sequence', 'annotation_method',
                     'exp_evidence', 'protein__name', 'protein__idtaxonomy__family',
                     'protein__idtaxonomy__commonname']
    pagination_class = CustomPagination


class HostAPIView(generics.ListCreateAPIView):
    queryset = Host.objects.all()
    serializer_class = HostSerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['idhost', 'host', 'ncbitax']
    search_fields = filterset_fields

    def get_object(self, pk=None):
        if pk is None: return
        try:
            return Host.objects.get(pk=pk)
        except Host.DoesNotExist:
            raise Http404

    def put(self, request, pk=None, format=None):
        if pk is None: return Response(status=status.HTTP_400_BAD_REQUEST)
        snippet = self.get_object(pk)
        serializer = HostSerializer(snippet, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk=None, format=None):
        if pk is None: return
        snippet = self.get_object(pk)
        snippet.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class InhibitorAntibodyAPIView(generics.ListCreateAPIView):
    queryset = InhibitorAntibody.objects.all()
    serializer_class = InhibitorAntibodySerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['idsubstance', 'type', 'repository', 'id_repository', 'idprotein']
    search_fields = ['idsubstance', 'type', 'repository', 'id_repository', 'idprotein__name',
                     'idprotein__idtaxonomy__commonname']

    def get_object(self, pk=None):
        if pk is None: return
        try:
            return InhibitorAntibody.objects.get(pk=pk)
        except InhibitorAntibody.DoesNotExist:
            raise Http404

    def put(self, request, pk=None, format=None):
        if pk is None: return Response(status=status.HTTP_400_BAD_REQUEST)
        snippet = self.get_object(pk)
        serializer = InhibitorAntibodySerializer(snippet, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk=None, format=None):
        if pk is None: return
        snippet = self.get_object(pk)
        snippet.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class InhibitorAntibodyAPIView_Save(generics.ListCreateAPIView):
    queryset = InhibitorAntibody.objects.all()
    serializer_class = InhibitorAntibodySerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['idsubstance', 'type', 'repository', 'id_repository', 'idprotein']
    search_fields = ['idsubstance', 'type', 'repository', 'id_repository', 'idprotein__name',
                     'idprotein__idtaxonomy__commonname']
    pagination_class = CustomPagination


class PeptideReferencesAPIView(generics.ListCreateAPIView):
    queryset = PeptideReferences.objects.all()
    serializer_class = PeptideReferencesSerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['idpeptide', 'idreferences', 'idreferences__doi',
                        'idpeptide__protein__idtaxonomy__commonname']
    search_fields = ['idpeptide__protein__name', 'idreferences__doi',
                     'idpeptide__protein__idtaxonomy__commonname']

    def get_object(self, pk=None):
        if pk is None: return
        try:
            return PeptideReferences.objects.get(pk=pk)
        except PeptideReferences.DoesNotExist:
            raise Http404

    def put(self, request, pk=None, format=None):
        if pk is None: return Response(status=status.HTTP_400_BAD_REQUEST)
        snippet = self.get_object(pk)
        serializer = PeptideReferencesSerializer(snippet, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk=None, format=None):
        if pk is None: return
        snippet = self.get_object(pk)
        snippet.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class PeptideReferencesAPIView_Save(generics.ListCreateAPIView):
    queryset = PeptideReferences.objects.all()
    serializer_class = PeptideReferencesSerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['idpeptide', 'idreferences', 'idreferences__doi',
                        'idpeptide__protein__idtaxonomy__commonname']
    search_fields = ['idpeptide__protein__name', 'idreferences__doi',
                     'idpeptide__protein__idtaxonomy__commonname']
    pagination_class = CustomPagination


class PeptideStructureAPIView(generics.ListCreateAPIView):
    queryset = PeptideStructure.objects.all()
    serializer_class = PeptideStructureSerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['idpeptide', 'idstructure',
                        'idstructure__id_repository']
    # search_fields = filterset_fields
    search_fields = ['idpeptide__protein__name', 'idpeptide__protein__idtaxonomy__commonname',
                     'idstructure__id_repository']

    def get_object(self, pk=None):
        if pk is None: return
        try:
            return PeptideStructure.objects.get(pk=pk)
        except PeptideStructure.DoesNotExist:
            raise Http404

    def put(self, request, pk=None, format=None):
        if pk is None: return Response(status=status.HTTP_400_BAD_REQUEST)
        snippet = self.get_object(pk)
        serializer = PeptideReferencesSerializer(snippet, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk=None, format=None):
        if pk is None: return
        snippet = self.get_object(pk)
        snippet.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class PeptideStructureAPIView_Save(generics.ListCreateAPIView):
    queryset = PeptideStructure.objects.all()
    serializer_class = PeptideStructureSerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['idpeptide', 'idstructure',
                        'idstructure__id_repository']
    # search_fields = filterset_fields
    search_fields = ['idpeptide__protein__name', 'idpeptide__protein__idtaxonomy__commonname',
                     'idstructure__id_repository']
    pagination_class = CustomPagination


class ProteinAPIView(generics.ListCreateAPIView):
    queryset = Protein.objects.all()
    serializer_class = ProteinSerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['idprotein', 'name', 'class_field', 'activation',
                        'name_fusogenic_unit', 'sequence_fusogenic',
                        'uniprotid', 'ncbiid', 'idtaxonomy', 'idtaxonomy__commonname', 'idtaxonomy__family']
    # search_fields = filterset_fields
    search_fields = ['idprotein', 'name', 'class_field', 'activation',
                     'name_fusogenic_unit', 'sequence_fusogenic',
                     'uniprotid', 'ncbiid', 'idtaxonomy__commonname', 'idtaxonomy__family']

    def get_object(self, pk=None):
        if pk is None: return
        try:
            return Protein.objects.get(pk=pk)
        except Protein.DoesNotExist:
            raise Http404

    def put(self, request, pk=None, format=None):
        if pk is None: return Response(status=status.HTTP_400_BAD_REQUEST)
        snippet = self.get_object(pk)
        serializer = ProteinSerializer(snippet, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk=None, format=None):
        if pk is None: return
        snippet = self.get_object(pk)
        snippet.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class ProteinAPIView_Save(generics.ListCreateAPIView):
    queryset = Protein.objects.all()
    serializer_class = ProteinSerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['idprotein', 'name', 'class_field', 'activation',
                        'name_fusogenic_unit', 'sequence_fusogenic',
                        'uniprotid', 'ncbiid', 'idtaxonomy', 'idtaxonomy__commonname']
    # search_fields = filterset_fields
    search_fields = ['idprotein', 'name', 'class_field', 'activation',
                     'name_fusogenic_unit', 'sequence_fusogenic',
                     'uniprotid', 'ncbiid', 'idtaxonomy__commonname']
    pagination_class = CustomPagination


class ProteinReferencesAPIView(generics.ListCreateAPIView):
    queryset = ProteinReferences.objects.all()
    serializer_class = ProteinReferencesSerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['idprotein', 'idreferences', 'idreferences__doi',
                        'idprotein__idtaxonomy__commonname']
    # search_fields = filterset_fields
    search_fields = ['idreferences__doi', 'idprotein__name',
                     'idprotein__idtaxonomy__commonname']

    def get_object(self, pk=None):
        if pk is None: return
        try:
            return ProteinReferences.objects.get(pk=pk)
        except ProteinReferences.DoesNotExist:
            raise Http404

    def put(self, request, pk=None, format=None):
        if pk is None: return Response(status=status.HTTP_400_BAD_REQUEST)
        snippet = self.get_object(pk)
        serializer = ProteinReferencesSerializer(snippet, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk=None, format=None):
        if pk is None: return
        snippet = self.get_object(pk)
        snippet.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class ProteinReferencesAPIView_Save(generics.ListCreateAPIView):
    queryset = ProteinReferences.objects.all()
    serializer_class = ProteinReferencesSerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['idprotein', 'idreferences', 'idreferences__doi',
                        'idprotein__idtaxonomy__commonname']
    # search_fields = filterset_fields
    search_fields = ['idreferences__doi', 'idprotein__name',
                     'idprotein__idtaxonomy__commonname']
    pagination_class = CustomPagination


class StructureAPIView(generics.ListCreateAPIView):
    queryset = Structure.objects.all()
    serializer_class = StructureSerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['idstructure', 'exp_method', 'repository', 'id_repository', 'reference__doi']
    search_fields = filterset_fields

    def get_object(self, pk=None):
        if pk is None: return
        try:
            return Structure.objects.get(pk=pk)
        except Structure.DoesNotExist:
            raise Http404

    def put(self, request, pk=None, format=None):
        if pk is None: return Response(status=status.HTTP_400_BAD_REQUEST)
        snippet = self.get_object(pk)
        serializer = StructureSerializer(snippet, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk=None, format=None):
        if pk is None: return
        snippet = self.get_object(pk)
        snippet.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class ReferencesAPIView(generics.ListCreateAPIView):
    queryset = References.objects.all()
    serializer_class = ReferencesSerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['idreferences', 'doi']
    search_fields = filterset_fields

    def get_object(self, pk=None):
        if pk is None: return
        try:
            return References.objects.get(pk=pk)
        except References.DoesNotExist:
            raise Http404

    def put(self, request, pk=None, format=None):
        if pk is None: return Response(status=status.HTTP_400_BAD_REQUEST)
        snippet = self.get_object(pk)
        serializer = ReferencesSerializer(snippet, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk=None, format=None):
        if pk is None: return
        snippet = self.get_object(pk)
        snippet.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class TaxHostAPIView(generics.ListCreateAPIView):
    queryset = TaxHost.objects.all()
    serializer_class = TaxHostSerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['idtaxonomy', 'idtaxonomy__commonname', 'idtaxonomy__ncbitax',
                        'idhost', 'idhost__host', 'idhost__ncbitax']
    # search_fields = filterset_fields
    search_fields = ['idtaxonomy__commonname', 'idtaxonomy__ncbitax',
                     'idhost__host', 'idhost__ncbitax']

    def get_object(self, pk=None):
        if pk is None: return
        try:
            return TaxHost.objects.get(pk=pk)
        except TaxHost.DoesNotExist:
            raise Http404

    def put(self, request, pk=None, format=None):
        if pk is None: return Response(status=status.HTTP_400_BAD_REQUEST)
        snippet = self.get_object(pk)
        serializer = TaxHostSerializer(snippet, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk=None, format=None):
        if pk is None: return
        snippet = self.get_object(pk)
        snippet.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class TaxHostAPIView_Save(generics.ListCreateAPIView):
    queryset = TaxHost.objects.all()
    serializer_class = TaxHostSerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['idtaxonomy', 'idtaxonomy__commonname', 'idtaxonomy__ncbitax',
                        'idhost', 'idhost__host', 'idhost__ncbitax']
    # search_fields = filterset_fields
    search_fields = ['idtaxonomy__commonname', 'idtaxonomy__ncbitax',
                     'idhost__host', 'idhost__ncbitax']
    pagination_class = CustomPagination


class TaxonomyVirusAPIView(generics.ListCreateAPIView):
    queryset = TaxonomyVirus.objects.all()
    serializer_class = TaxonomyVirusSerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['idtaxonomy', 'commonname', 'family', 'genre',
                        'species', 'subspecies', 'ncbitax']
    search_fields = filterset_fields

    def get_object(self, pk=None):
        if pk is None: return
        try:
            return TaxonomyVirus.objects.get(pk=pk)
        except TaxonomyVirus.DoesNotExist:
            raise Http404

    def put(self, request, pk=None, format=None):
        if pk is None: return Response(status=status.HTTP_400_BAD_REQUEST)
        snippet = self.get_object(pk)
        serializer = TaxonomyVirusSerializer(snippet, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk=None, format=None):
        if pk is None: return
        snippet = self.get_object(pk)
        snippet.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class TaxonomyVirusAPIView_Save(generics.ListCreateAPIView):
    queryset = TaxonomyVirus.objects.all()
    serializer_class = TaxonomyVirusSerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['idtaxonomy', 'commonname', 'family', 'genre',
                        'species', 'subspecies', 'ncbitax']
    search_fields = filterset_fields
    pagination_class = CustomPagination
