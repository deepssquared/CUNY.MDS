CREATE TABLE participants
(first_name CHAR(25),
 last_name CHAR(25),
 participant_id INTEGER(1),
 )

CREATE INDEX part_id
ON participants (participant_id);

INSERT INTO participants (first_name, last_name, participant_id)
VALUES ('PEDRO', 'BIRDSALL', 1);

INSERT INTO participants (first_name, last_name, participant_id)
VALUES ('JOE', 'EXOCTIC', 2);

INSERT INTO participants (first_name, last_name, participant_id)
VALUES ('CAROLE', 'BASKIN', 3);

INSERT INTO participants (first_name, last_name, participant_id)
VALUES ('MEGHAN', 'MARKLE', 4);

INSERT INTO participants (first_name, last_name, participant_id)
VALUES ('LOVE', 'QUINN', 5);

CREATE TABLE ratings
(participant_id INTEGER(1),
movie CHAR(50),
rating INTEGER(1),
);
 
CREATE INDEX part_id
ON ratings (participant_id);

INSERT INTO ratings (participant_id, movie, rating)
VALUES (1, 'SPIDERMAN', 2);

INSERT INTO ratings (participant_id, movie, rating)
VALUES (1, 'SCREAM', 2);

INSERT INTO ratings (participant_id, movie, rating)
VALUES (1, 'WEST SIDE STORY', 4);

INSERT INTO ratings (participant_id, movie, rating)
VALUES (2, 'ENCANTO', 5);

INSERT INTO ratings (participant_id, movie, rating)
VALUES (2, 'WEST SIDE STORY', 2);

INSERT INTO ratings (participant_id, movie, rating)
VALUES (3, 'THE LAST DUEL', 1);

INSERT INTO ratings (participant_id, movie, rating)
VALUES (3, 'SCREAM', 5);

INSERT INTO ratings (participant_id, movie, rating)
VALUES (3, 'ENCANTO', 1);

INSERT INTO ratings (participant_id, movie, rating)
VALUES (4, 'ENCANTO', 3);

INSERT INTO ratings (participant_id, movie, rating)
VALUES (4, 'THE LAST DUEL', 5);

INSERT INTO ratings (participant_id, movie, rating)
VALUES (5, 'LICORICE PIZZA', 3);
