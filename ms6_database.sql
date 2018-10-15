DROP DATABASE ms6;

CREATE DATABASE ms6;

USE ms6;

/** CRIA A TABELA DE REQUISICOES ("jobs") **/

CREATE TABLE jobs (id MEDIUMINT NOT NULL AUTO_INCREMENT PRIMARY KEY,
                   job_id CHAR(100),
                   email CHAR(100),
                   status CHAR(100),
                   replicates CHAR(100),
                   reference_type CHAR(10),
                   idioma CHAR(10),
                   peptideshaker_parameters CHAR(100),
                   contaminants CHAR(100)
                   );