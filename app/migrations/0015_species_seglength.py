# Generated by Django 4.1.3 on 2023-12-14 02:48

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('app', '0014_rename_sample_id_position_sample'),
    ]

    operations = [
        migrations.AddField(
            model_name='species',
            name='seglength',
            field=models.TextField(null=True),
        ),
    ]
