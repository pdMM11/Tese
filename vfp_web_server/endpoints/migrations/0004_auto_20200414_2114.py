# Generated by Django 3.0.3 on 2020-04-14 20:14

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('endpoints', '0003_auto_20200414_2112'),
    ]

    operations = [
        migrations.AlterField(
            model_name='abtest',
            name='summary',
            field=models.TextField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='abtest',
            name='title',
            field=models.TextField(),
        ),
    ]
