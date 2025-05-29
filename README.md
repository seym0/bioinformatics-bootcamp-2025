# Анализ геномных данных дрожжей после модификации

## Описание 
Проект был предложен в рамках [хакатона по прикладной геномике Bioinformatics Bootcamp 2025 от ПИШ ИТМО](https://pish.itmo.ru/genomics-bootcamp). Проект включает обработку данных WGS дрожжей, подвергнутых CRISPR-редактированию. Результаты важны для биотехнологии и синтетической биологии.

## Данные

- **Референсный геном**: *Saccharomyces cerevisiae* S288C (`GCF_000146045.2_R64_genomic.fna`).
- **Сэмплы**: 16 CRISPR-редактированных штаммов (`reads/Sample_hst3.4.1` – `reads/Sample_hst3.4.16`) (мутанты).
- **Контроль**: Немодифицированный штамм *S288C-WT* (`reads/Sample_control`).
- **Формат данных**: FASTQ-файлы (`*.fastq.gz`)

# Решение
## Проверка качества сырых прочтений

Для проверки качества прочтений семплов данных парного секвенирования с двумя лейнами мы использовали [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) в пайплайне ```Snakefile_fastqc```, а затем получившиеся отчеты были объединены в общий отчет ```multiqc_report_raw_reads.html``` с помощью [multiqc](https://seqera.io/multiqc/).

## Обрезка сырых прочтений

После ручного анализа качества сырых прочтений для удаления Overrepresented sequences в первом, втором, седьмом и двенадцатом семлпе использовали [fastp](https://github.com/OpenGene/fastp) в пайплайне ```Snakefile_fastp```. Обработанные сэмплы и исходные нетронутые объединили в новую директорию ```preprocessed_reads/```.

## Выравнивание на референсный геном

Для индексации [референсного генома S.cerevisiae S288C](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/) использовали [bwa](https://github.com/lh3/bwa): 

```bash
bwa index S288C/GCF_000146045.2_R64_genomic.fna
```

Затем выполнили выравнивание на референсный геном для всех 16 семплов и контроля - генома дрожжей . Для каждого семпла обработали две пары ридов (Lane 1: R1+R2, Lane 2: R1+R2; лейн — это техническое разделение данных одного семпла), выровняли на референс и объединили результаты выравнивания в один BAM-файл на семпл и отсортировали с помощью [samtools](https://www.htslib.org/). 

Note: Объединение BAM-файлов (а не конкатенация fastq файлов для прямого и обратного прочтения из лейнов перед выравнивания) — более стандартный подход, так как позволяет сохранить метаданные (например, Read Groups) для каждого лейна. 

```bash
mkdir -p aligned

for i in {1..16}; do
    echo "Processing sample${i}..."
    # Lane 1
    l1_r1=$(ls ./preprocessed_reads/Sample_hst3.4.${i}/hst3.4.${i}_*L001_R1*.fastq.gz)
    l1_r2=$(ls ./preprocessed_reads/Sample_hst3.4.${i}/hst3.4.${i}_*L001_R2*.fastq.gz)
    
    bwa mem -t 4 -R "@RG\tID:sample${i}_L001\tSM:sample${i}\tLB:lib1\tPL:ILLUMINA" ../Samples/S288C/GCF_000146045.2_R64_genomic.fna "$l1_r1" "$l1_r2" | samtools view -bS - | samtools sort -o aligned/sample${i}_L001.bam

    # Lane 2
    l2_r1=$(ls ./preprocessed_reads/Sample_hst3.4.${i}/hst3.4.${i}_*L002_R1*.fastq.gz)
    l2_r2=$(ls ./preprocessed_reads/Sample_hst3.4.${i}/hst3.4.${i}_*L002_R2*.fastq.gz)
    
    bwa mem -t 4 -R "@RG\tID:sample${i}_L002\tSM:sample${i}\tLB:lib1\tPL:ILLUMINA" ../Samples/S288C/GCF_000146045.2_R64_genomic.fna "$l2_r1" "$l2_r2" | samtools view -bS - | samtools sort -o aligned/sample${i}_L002.bam
    
    # Объединение BAM-файлов и их индексация 
    samtools merge -@ 8 aligned/sorted_sample${i}.bam aligned/sample${i}_L001.bam aligned/sample${i}_L002.bam
    samtools index aligned/sorted_sample${i}.bam

    # Можно удалить файлы после выравнивания и их бинарники, если потребуется
    # rm aligned/sample${i}_L001.bam aligned/sample${i}_L002.bam aligned/sample${i}_*.sam
done
```

## Генотипирование


## Аннотация