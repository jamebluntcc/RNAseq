'''
Usage:
    add_sp_inf.py <annotation.cfg>


'''
from sqlalchemy import create_engine, Column, Integer, String, Sequence, DateTime
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import and_
import datetime
from ConfigParser import ConfigParser
from docopt import docopt

Base = declarative_base()


class species_annotation_info(Base):
    __tablename__ = 'species_annotation_info'

    id = Column(Integer, Sequence('user_id_seq'),
                primary_key=True, autoincrement=True)

    species_latin = Column(String(100))
    species_database = Column(String(50))
    species_database_version = Column(String(50))
    upload_date = Column(DateTime, default=datetime.datetime.utcnow)

    def __repr__(self):
        return "{0}|{1}|{2}|{3}".format(self.species_latin, self.species_database, self.species_database_version, self.upload_date)

if __name__ == '__main__':
    eng = create_engine('mysql://lxgui:lxgui@localhost/test1')
    arguments = docopt(__doc__, version='v1')
    ann_cfg = arguments['<annotation.cfg>']
    conf = ConfigParser()
    conf.read(ann_cfg)
    species_latin = conf.get('annotation', 'species_latin')
    species_database = conf.get('annotation', 'species_database')
    species_database_version = conf.get(
        'annotation', 'species_database_version')

    if not eng.dialect.has_table(eng, 'species_annotation_info'):
        Base.metadata.create_all(bind=eng)
    Session = sessionmaker(bind=eng)
    session = Session()
    rs = session.query(species_annotation_info).filter(and_(species_annotation_info.species_latin == species_latin,
                                                            species_annotation_info.species_database == species_database,
                                                            species_annotation_info.species_database_version == species_database_version))
    if not rs.count():
        session.add(species_annotation_info(species_latin=species_latin,
                                            species_database=species_database,
                                            species_database_version=species_database_version))
    session.commit()
