
/*
RABDAM
Copyright (C) 2018 Garman Group, University of Oxford

This file is part of RABDAM.

RABDAM is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version.

RABDAM is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General
Public License along with this program.  If not, see
<http://www.gnu.org/licenses/>.
*/

$(document).ready(function() {
    /* Hides tables from view */
    $('table').hide();

    /* Underlines sublist headings when mouse is moved over them */
    $('h3').hover(function() {
        $(this).css('text-decoration', 'underline');
     }, function(){
        $(this).css('text-decoration', 'none');
      });

    /* Changes mouse from cursor to pointer when mouse is moved over sublist
     headings */
     $('h3').hover(function() {
         $(this).css('cursor', 'pointer');
      }, function(){
         $(this).css('cursor', 'default');
       });

    /* Reveals tables when relevant sublist heading is clicked */
    $('.sublist').click(function() {
        $(this).find('table').slideToggle('fast');
    });

});
